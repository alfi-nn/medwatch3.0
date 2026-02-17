"""Improved ADR prediction with focal loss, fewer classes, and better training."""
import torch
import torch.nn as nn
import torch.nn.functional as F
from pathlib import Path
from torch.optim import AdamW
from torch.optim.lr_scheduler import CosineAnnealingLR
from sklearn.metrics import roc_auc_score
import numpy as np
from torch_geometric.nn import HGTConv, Linear

DRIVE_PROCESSED = Path('/content/drive/MyDrive/MedWatch3/data/processed')
DRIVE_CHECKPOINTS = Path('/content/drive/MyDrive/MedWatch3/checkpoints')
device = 'cuda' if torch.cuda.is_available() else 'cpu'

# === Use Top 20 Side Effects Only (easier task) ===
NUM_SIDE_EFFECTS = 20
EPOCHS = 300
LEARNING_RATE = 0.001

class FocalLoss(nn.Module):
    """Focal loss for imbalanced classification."""
    def __init__(self, alpha=0.25, gamma=2.0):
        super().__init__()
        self.alpha = alpha
        self.gamma = gamma
    
    def forward(self, pred, target):
        bce = F.binary_cross_entropy(pred, target, reduction='none')
        pt = torch.where(target == 1, pred, 1 - pred)
        focal_weight = self.alpha * (1 - pt) ** self.gamma
        return (focal_weight * bce).mean()

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
        self.fc1 = Linear(hidden_channels, hidden_channels * 2)
        self.fc2 = Linear(hidden_channels * 2, hidden_channels)
        self.fc3 = Linear(hidden_channels, num_side_effects)
        self.dropout = nn.Dropout(0.3)
    
    def forward(self, x):
        x = F.relu(self.fc1(x))
        x = self.dropout(x)
        x = F.relu(self.fc2(x))
        x = self.dropout(x)
        return torch.sigmoid(self.fc3(x))

def train():
    # Load data
    data = torch.load(DRIVE_PROCESSED / 'graph_data.pt', weights_only=False)
    labels_data = torch.load(DRIVE_PROCESSED / 'side_effect_labels.pt', weights_only=False)
    
    # Use only top N side effects
    labels = labels_data['labels'][:, :NUM_SIDE_EFFECTS].to(device)
    side_effects = labels_data['side_effects'][:NUM_SIDE_EFFECTS]
    
    print(f"Side effects: {NUM_SIDE_EFFECTS}")
    print(f"Side effect names: {side_effects}")
    print(f"Labels shape: {labels.shape}")
    print(f"Positive labels per SE: {labels.sum(dim=0).tolist()}")
    
    in_channels = {nt: data[nt].x.shape[1] for nt in data.node_types}
    
    # Models - train HGT from scratch for ADR task
    hgt = HGT(in_channels, 256, 4, 2, data.metadata()).to(device)
    adr_pred = ADRPredictor(256, NUM_SIDE_EFFECTS).to(device)
    
    # Focal loss for imbalanced classes
    criterion = FocalLoss(alpha=0.25, gamma=2.0)
    
    # AdamW with weight decay
    optimizer = AdamW(
        list(hgt.parameters()) + list(adr_pred.parameters()), 
        lr=LEARNING_RATE, 
        weight_decay=0.01
    )
    
    # Cosine annealing scheduler
    scheduler = CosineAnnealingLR(optimizer, T_max=EPOCHS)
    
    x_dict = {k: v.to(device) for k, v in data.x_dict.items()}
    edge_index_dict = {k: v.to(device) for k, v in data.edge_index_dict.items()}
    
    # Train/val split
    n_drugs = labels.shape[0]
    perm = torch.randperm(n_drugs)
    train_idx = perm[:int(0.8 * n_drugs)]
    val_idx = perm[int(0.8 * n_drugs):]
    
    best_auc = 0
    patience = 50
    patience_counter = 0
    
    for epoch in range(1, EPOCHS + 1):
        hgt.train()
        adr_pred.train()
        
        z_dict = hgt(x_dict, edge_index_dict)
        drug_emb = z_dict['drug']
        
        pred = adr_pred(drug_emb[train_idx])
        loss = criterion(pred, labels[train_idx])
        
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
                pred_val = adr_pred(z_dict['drug'][val_idx])
                
                pred_np = pred_val.cpu().numpy()
                labels_np = labels[val_idx].cpu().numpy()
                
                aucs = []
                for i in range(NUM_SIDE_EFFECTS):
                    if labels_np[:, i].sum() > 0 and labels_np[:, i].sum() < len(labels_np):
                        try:
                            auc = roc_auc_score(labels_np[:, i], pred_np[:, i])
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
                }, DRIVE_CHECKPOINTS / 'best_adr_improved.pt')
                print(f"  💾 Saved! (AUC: {avg_auc:.4f})")
            else:
                patience_counter += 10
                if patience_counter >= patience:
                    print(f"Early stopping at epoch {epoch}")
                    break
    
    print(f"\n✅ Training complete! Best AUC: {best_auc:.4f}")
    print(f"Model saved to: {DRIVE_CHECKPOINTS / 'best_adr_improved.pt'}")

if __name__ == "__main__":
    train()
