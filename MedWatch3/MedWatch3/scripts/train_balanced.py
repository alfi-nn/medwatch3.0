"""
Phase 3: Train ADR model on balanced/augmented data.
Architecture: HGT(256, 4 heads, 2 layers) + ADRPredictor(256->512->256->N)
Loss: BCEWithLogitsLoss with per-class pos_weight
"""
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

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / 'data' / 'processed'
CHECKPOINT_DIR = BASE_DIR / 'checkpoints'
DEVICE = 'cuda' if torch.cuda.is_available() else 'cpu'

EPOCHS = 200
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
    """Same 3-layer architecture as best_adr_improved."""
    def __init__(self, hidden_channels, num_side_effects):
        super().__init__()
        self.fc1 = Linear(hidden_channels, hidden_channels * 2)   # 256 -> 512
        self.fc2 = Linear(hidden_channels * 2, hidden_channels)    # 512 -> 256
        self.fc3 = Linear(hidden_channels, num_side_effects)       # 256 -> N
        self.dropout = nn.Dropout(0.3)
    
    def forward(self, x):
        x = F.relu(self.fc1(x))
        x = self.dropout(x)
        x = F.relu(self.fc2(x))
        x = self.dropout(x)
        return self.fc3(x)  # Raw logits (no sigmoid, BCEWithLogitsLoss handles it)

def train():
    CHECKPOINT_DIR.mkdir(exist_ok=True)
    
    print(f"Training on {DEVICE}...")
    
    # Load augmented data
    data = torch.load(DATA_DIR / 'graph_data_augmented.pt', weights_only=False, map_location='cpu')
    labels_data = torch.load(DATA_DIR / 'side_effect_labels_augmented.pt', weights_only=False, map_location='cpu')
    
    labels = labels_data['labels'].to(DEVICE)
    side_effects = labels_data['side_effects']
    n_original = labels_data['n_original']
    num_se = len(side_effects)
    
    print(f"Drugs: {labels.shape[0]} ({n_original} original + {labels_data['n_augmented']} augmented)")
    print(f"Side effects: {num_se}")
    print(f"Mean positive rate: {labels.mean():.2%}")
    
    # Calculate pos_weight for class imbalance
    pos_counts = labels.sum(dim=0)
    neg_counts = labels.shape[0] - pos_counts
    pos_weights = (neg_counts / (pos_counts + 1e-6)).to(DEVICE)
    print(f"Avg pos_weight: {pos_weights.mean():.2f}")
    
    in_channels = {nt: data[nt].x.shape[1] for nt in data.node_types}
    
    # Models  
    hgt = HGT(in_channels, HIDDEN_CHANNELS, 4, 2, data.metadata()).to(DEVICE)
    adr_pred = ADRPredictor(HIDDEN_CHANNELS, num_se).to(DEVICE)
    
    # Loss with pos_weight
    criterion = nn.BCEWithLogitsLoss(pos_weight=pos_weights)
    
    optimizer = AdamW(
        list(hgt.parameters()) + list(adr_pred.parameters()), 
        lr=LEARNING_RATE, 
        weight_decay=0.01
    )
    scheduler = CosineAnnealingLR(optimizer, T_max=EPOCHS)
    
    x_dict = {k: v.to(DEVICE) for k, v in data.x_dict.items()}
    edge_index_dict = {k: v.to(DEVICE) for k, v in data.edge_index_dict.items()}
    
    # Split: use original drugs for val, augmented for train
    # This prevents data leakage from augmented copies
    n_total = labels.shape[0]
    
    # Use 80% of original drugs for train, 20% for val
    # All augmented drugs go to train
    perm = torch.randperm(n_original)
    n_val = int(0.2 * n_original)
    val_original_idx = perm[:n_val]
    train_original_idx = perm[n_val:]
    
    # Augmented indices start at n_original
    augmented_idx = torch.arange(n_original, n_total)
    
    train_idx = torch.cat([train_original_idx, augmented_idx])
    val_idx = val_original_idx
    
    print(f"Train: {len(train_idx)} drugs, Val: {len(val_idx)} drugs")
    
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
                    col = labels_np[:, i]
                    if 0 < col.sum() < len(col):
                        try:
                            auc = roc_auc_score(col, probs_np[:, i])
                            aucs.append(auc)
                        except:
                            pass
                
                avg_auc = np.mean(aucs) if aucs else 0
            
            lr = scheduler.get_last_lr()[0]
            print(f"Epoch {epoch:03d} | Loss: {loss:.4f} | Val AUC: {avg_auc:.4f} | LR: {lr:.6f}")
            
            if avg_auc > best_auc:
                best_auc = avg_auc
                patience_counter = 0
                torch.save({
                    'epoch': epoch,
                    'hgt_state_dict': hgt.state_dict(),
                    'adr_state_dict': adr_pred.state_dict(),
                    'auc': avg_auc,
                    'side_effects': side_effects,
                    'num_layers': 3,
                    'hidden_channels': HIDDEN_CHANNELS,
                }, CHECKPOINT_DIR / 'best_adr_balanced_v2.pt')
                print(f"  Saved! (AUC: {avg_auc:.4f})")
            else:
                patience_counter += 1
                if patience_counter >= patience:
                    print(f"Early stopping at epoch {epoch}")
                    break
    
    print(f"\nTraining complete! Best AUC: {best_auc:.4f}")
    print(f"Model saved to: {CHECKPOINT_DIR / 'best_adr_balanced_v2.pt'}")

if __name__ == "__main__":
    train()
