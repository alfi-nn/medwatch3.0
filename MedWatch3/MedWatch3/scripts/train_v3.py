"""
V3 Training: Optimized from scratch with label smoothing, warm-up, and stronger regularization.
"""
import torch
import torch.nn as nn
import torch.nn.functional as F
from pathlib import Path
from torch.optim import AdamW
from sklearn.metrics import roc_auc_score, average_precision_score
import numpy as np
from torch_geometric.nn import HGTConv, Linear

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / 'data' / 'processed'
CHECKPOINT_DIR = BASE_DIR / 'checkpoints'
DEVICE = 'cpu'

EPOCHS = 300
LEARNING_RATE = 0.002
HIDDEN_CHANNELS = 256
LABEL_SMOOTH = 0.05  # Smoothing factor
WARMUP_EPOCHS = 20

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
        self.bn1 = nn.BatchNorm1d(hidden_channels * 2)
        self.fc2 = Linear(hidden_channels * 2, hidden_channels)
        self.bn2 = nn.BatchNorm1d(hidden_channels)
        self.fc3 = Linear(hidden_channels, num_side_effects)
        self.dropout = nn.Dropout(0.4)
    
    def forward(self, x):
        x = F.relu(self.bn1(self.fc1(x)))
        x = self.dropout(x)
        x = F.relu(self.bn2(self.fc2(x)))
        x = self.dropout(x)
        return self.fc3(x)

def get_lr(epoch, warmup, base_lr, total_epochs):
    """Warm-up + cosine decay."""
    if epoch < warmup:
        return base_lr * (epoch + 1) / warmup
    progress = (epoch - warmup) / (total_epochs - warmup)
    return base_lr * 0.5 * (1 + np.cos(np.pi * progress))

def train():
    CHECKPOINT_DIR.mkdir(exist_ok=True)
    print("V3 Training: Label Smoothing + BatchNorm + Warm-up")
    
    # Load augmented data
    data = torch.load(DATA_DIR / 'graph_data_augmented.pt', weights_only=False, map_location='cpu')
    labels_data = torch.load(DATA_DIR / 'side_effect_labels_augmented.pt', weights_only=False, map_location='cpu')
    
    labels = labels_data['labels'].to(DEVICE)
    side_effects = labels_data['side_effects']
    n_original = labels_data['n_original']
    num_se = len(side_effects)
    
    # Label smoothing: convert 0/1 to 0.05/0.95
    labels_smooth = labels * (1 - 2 * LABEL_SMOOTH) + LABEL_SMOOTH
    
    print(f"Drugs: {labels.shape[0]} ({n_original} original + {labels_data['n_augmented']} augmented)")
    print(f"Side effects: {num_se}")
    
    # Calculate pos_weight
    pos_counts = labels.sum(dim=0)
    neg_counts = labels.shape[0] - pos_counts
    pos_weights = (neg_counts / (pos_counts + 1e-6)).clamp(0.5, 2.0).to(DEVICE)
    
    in_channels = {nt: data[nt].x.shape[1] for nt in data.node_types}
    
    hgt = HGT(in_channels, HIDDEN_CHANNELS, 4, 2, data.metadata()).to(DEVICE)
    adr_pred = ADRPredictor(HIDDEN_CHANNELS, num_se).to(DEVICE)
    
    criterion = nn.BCEWithLogitsLoss(pos_weight=pos_weights)
    
    all_params = list(hgt.parameters()) + list(adr_pred.parameters())
    optimizer = AdamW(all_params, lr=LEARNING_RATE, weight_decay=0.01)
    
    x_dict = {k: v.to(DEVICE) for k, v in data.x_dict.items()}
    edge_index_dict = {k: v.to(DEVICE) for k, v in data.edge_index_dict.items()}
    
    # Split
    perm = torch.randperm(n_original)
    n_val = int(0.2 * n_original)
    val_idx = perm[:n_val]
    train_original_idx = perm[n_val:]
    augmented_idx = torch.arange(n_original, labels.shape[0])
    train_idx = torch.cat([train_original_idx, augmented_idx])
    
    print(f"Train: {len(train_idx)}, Val: {len(val_idx)}")
    
    best_auc = 0
    patience = 40
    patience_counter = 0
    
    for epoch in range(1, EPOCHS + 1):
        # Warm-up + cosine LR
        lr = get_lr(epoch, WARMUP_EPOCHS, LEARNING_RATE, EPOCHS)
        for pg in optimizer.param_groups:
            pg['lr'] = lr
        
        hgt.train()
        adr_pred.train()
        
        z_dict = hgt(x_dict, edge_index_dict)
        logits = adr_pred(z_dict['drug'][train_idx])
        loss = criterion(logits, labels_smooth[train_idx])
        
        optimizer.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(all_params, 1.0)
        optimizer.step()
        
        if epoch % 10 == 0:
            hgt.eval()
            adr_pred.eval()
            with torch.no_grad():
                z_dict = hgt(x_dict, edge_index_dict)
                logits_val = adr_pred(z_dict['drug'][val_idx])
                probs_val = torch.sigmoid(logits_val)
                
                probs_np = probs_val.cpu().numpy()
                labels_np = labels[val_idx].cpu().numpy()  # Use unsmoothed for eval
                
                aucs = []
                for i in range(num_se):
                    col = labels_np[:, i]
                    if 0 < col.sum() < len(col):
                        try:
                            aucs.append(roc_auc_score(col, probs_np[:, i]))
                        except:
                            pass
                
                avg_auc = np.mean(aucs) if aucs else 0
            
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
                    'has_batchnorm': True,
                }, CHECKPOINT_DIR / 'best_adr_v3.pt')
                print(f"  Saved! (AUC: {avg_auc:.4f})")
            else:
                patience_counter += 1
                if patience_counter >= patience:
                    print(f"Early stopping at epoch {epoch}")
                    break
    
    print(f"\nTraining complete! Best AUC: {best_auc:.4f}")

if __name__ == "__main__":
    train()
