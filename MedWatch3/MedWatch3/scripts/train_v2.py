
import torch
import torch.nn as nn
import torch.nn.functional as F
from pathlib import Path
from torch.optim import AdamW
from torch.optim.lr_scheduler import CosineAnnealingLR
from sklearn.metrics import roc_auc_score
import numpy as np
from torch_geometric.nn import HGTConv, Linear
import sys

# Add project root to path
sys.path.append(str(Path(__file__).resolve().parent.parent))

# Configuration
DRIVE_PROCESSED = Path('data/processed')
DRIVE_CHECKPOINTS = Path('checkpoints')
DEVICE = 'cuda' if torch.cuda.is_available() else 'cpu'

# Use Top 100 Side Effects (more challenging but better coverage)
NUM_SIDE_EFFECTS = 100
EPOCHS = 300
LEARNING_RATE = 0.001
HIDDEN_CHANNELS = 256

class HGT(nn.Module):
    def __init__(self, in_channels_dict, hidden_channels, num_heads, num_layers, metadata):
        super().__init__()
        self.lin_dict = nn.ModuleDict()
        for node_type in metadata[0]:
            self.lin_dict[node_type] = Linear(in_channels_dict[node_type], hidden_channels)
        self.convs = nn.ModuleList()
        for _ in range(num_layers):
            self.convs.append(HGTConv(hidden_channels, hidden_channels, metadata, num_heads))
    
    def forward(self, x_dict, edge_index_dict):
        out = {nt: F.dropout(self.lin_dict[nt](x).relu(), p=0.1, training=self.training) 
               for nt, x in x_dict.items()}
        for conv in self.convs:
            out = conv(out, edge_index_dict)
            out = {k: F.dropout(v, p=0.1, training=self.training) for k, v in out.items()}
        return out

class ADRPredictor(nn.Module):
    def __init__(self, hidden_channels, num_side_effects):
        super().__init__()
        # 3-layer MLP for better capacity
        self.net = nn.Sequential(
            Linear(hidden_channels, 512), nn.ReLU(), nn.Dropout(0.3),
            Linear(512, 256), nn.ReLU(), nn.Dropout(0.3),
            Linear(256, num_side_effects) # No sigmoid here, using BCEWithLogitsLoss
        )
    
    def forward(self, x):
        return self.net(x)

def train():
    print(f"Training on {DEVICE}...")
    
    # Ensure directories exist
    DRIVE_CHECKPOINTS.mkdir(exist_ok=True)
    
    # Load data
    try:
        data = torch.load(DRIVE_PROCESSED / 'graph_data.pt', weights_only=False)
        labels_data = torch.load(DRIVE_PROCESSED / 'side_effect_labels.pt', weights_only=False)
    except FileNotFoundError:
        print("Error: processed data not found. Please run data processing scripts first.")
        return

    # Use top N side effects
    # Check if we have enough side effects
    total_se = labels_data['labels'].shape[1]
    num_se = min(NUM_SIDE_EFFECTS, total_se)
    
    labels = labels_data['labels'][:, :num_se].to(DEVICE)
    side_effects = labels_data['side_effects'][:num_se]
    
    print(f"Side effects: {num_se}")
    print(f"Labels shape: {labels.shape}")
    
    # Calculate positive weights for class imbalance
    # pos_weight = (num_neg) / (num_pos)
    num_pos = labels.sum(dim=0)
    num_neg = labels.shape[0] - num_pos
    pos_weights = (num_neg / (num_pos + 1e-6)).to(DEVICE) # Add epsilon to avoid div by zero
    
    print(f"Avg positive weight: {pos_weights.mean().item():.2f}")
    
    in_channels = {nt: data[nt].x.shape[1] for nt in data.node_types}
    
    # Models
    hgt = HGT(in_channels, HIDDEN_CHANNELS, 4, 2, data.metadata()).to(DEVICE)
    adr_pred = ADRPredictor(HIDDEN_CHANNELS, num_se).to(DEVICE)
    
    # Loss function with weights
    criterion = nn.BCEWithLogitsLoss(pos_weight=pos_weights)
    
    optimizer = AdamW(
        list(hgt.parameters()) + list(adr_pred.parameters()), 
        lr=LEARNING_RATE, 
        weight_decay=0.01
    )
    
    scheduler = CosineAnnealingLR(optimizer, T_max=EPOCHS)
    
    x_dict = {k: v.to(DEVICE) for k, v in data.x_dict.items()}
    edge_index_dict = {k: v.to(DEVICE) for k, v in data.edge_index_dict.items()}
    
    # Train/val split
    n_drugs = labels.shape[0]
    perm = torch.randperm(n_drugs)
    train_idx = perm[:int(0.8 * n_drugs)]
    val_idx = perm[int(0.8 * n_drugs):]
    
    best_auc = 0
    patience = 30
    patience_counter = 0
    
    for epoch in range(1, EPOCHS + 1):
        hgt.train()
        adr_pred.train()
        
        z_dict = hgt(x_dict, edge_index_dict)
        drug_emb = z_dict['drug']
        
        logits = adr_pred(drug_emb[train_idx])
        loss = criterion(logits, labels[train_idx])
        
        optimizer.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(list(hgt.parameters()) + list(adr_pred.parameters()), 1.0)
        optimizer.step()
        scheduler.step()
        
        if epoch % 10 == 0:
            hgt.eval()
            adr_pred.eval()
            with torch.no_grad():
                z_dict = hgt(x_dict, edge_index_dict)
                logits_val = adr_pred(z_dict['drug'][val_idx])
                probs_val = torch.sigmoid(logits_val)
                
                probs_np = probs_val.cpu().numpy()
                labels_np = labels[val_idx].cpu().numpy()
                
                aucs = []
                for i in range(num_se):
                    if 0 < labels_np[:, i].sum() < len(labels_np):
                        try:
                            auc = roc_auc_score(labels_np[:, i], probs_np[:, i])
                            aucs.append(auc)
                        except:
                            pass
                
                avg_auc = np.mean(aucs) if aucs else 0
            
            lr = scheduler.get_last_lr()[0]
            print(f"Epoch {epoch:03d} | Loss: {loss:.4f} | Val AUC: {avg_auc:.4f} | LR: {lr:.6f}")
            
            if avg_auc > best_auc:
                best_auc = avg_auc
                patience_counter = 0
                
                # Save model
                state_dict = {
                    'epoch': epoch,
                    'hgt_state_dict': hgt.state_dict(),
                    'adr_state_dict': adr_pred.state_dict(), # We'll need to handle the 'net' key in load
                    'auc': avg_auc,
                    'side_effects': side_effects,
                }
                torch.save(state_dict, DRIVE_CHECKPOINTS / 'best_adr_v2.pt')
                print(f"  💾 Saved! (AUC: {avg_auc:.4f})")
            else:
                patience_counter += 1
                if patience_counter >= patience:
                    print(f"Early stopping at epoch {epoch}")
                    break
    
    print(f"\n✅ Training complete! Best AUC: {best_auc:.4f}")
    print(f"Model saved to: {DRIVE_CHECKPOINTS / 'best_adr_v2.pt'}")

if __name__ == "__main__":
    train()
