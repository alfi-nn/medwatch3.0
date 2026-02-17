"""Train on BALANCED side effects (50% positive rate)."""
import torch
import torch.nn as nn
import torch.nn.functional as F
from pathlib import Path
from torch.optim import AdamW
from sklearn.metrics import roc_auc_score
import numpy as np
from torch_geometric.nn import HGTConv, Linear

DRIVE_PROCESSED = Path('/content/drive/MyDrive/MedWatch3/data/processed')
DRIVE_CHECKPOINTS = Path('/content/drive/MyDrive/MedWatch3/checkpoints')
device = 'cuda' if torch.cuda.is_available() else 'cpu'

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
        out = {nt: F.dropout(F.relu(self.lin_dict[nt](x)), 0.1, self.training) 
               for nt, x in x_dict.items()}
        for conv in self.convs:
            out = conv(out, edge_index_dict)
        return out

class ADRPredictor(nn.Module):
    def __init__(self, hidden, num_se):
        super().__init__()
        self.net = nn.Sequential(
            Linear(hidden, hidden * 2), nn.ReLU(), nn.Dropout(0.3),
            Linear(hidden * 2, hidden), nn.ReLU(), nn.Dropout(0.3),
            Linear(hidden, num_se), nn.Sigmoid()
        )
    def forward(self, x):
        return self.net(x)

def train():
    data = torch.load(DRIVE_PROCESSED / 'graph_data.pt', weights_only=False)
    labels_data = torch.load(DRIVE_PROCESSED / 'side_effect_labels.pt', weights_only=False)
    
    all_labels = labels_data['labels']
    all_se = labels_data['side_effects']
    
    # Find BALANCED side effects (20-80% positive rate)
    n_drugs = all_labels.shape[0]
    pos_rates = all_labels.sum(dim=0) / n_drugs
    
    balanced_idx = []
    for i, rate in enumerate(pos_rates):
        if 0.2 <= rate <= 0.8:  # Between 20% and 80%
            balanced_idx.append(i)
    
    if len(balanced_idx) < 10:
        # Fallback: take side effects with most variance
        variances = all_labels.float().var(dim=0)
        balanced_idx = variances.argsort(descending=True)[:20].tolist()
    
    balanced_idx = balanced_idx[:20]  # Max 20
    
    labels = all_labels[:, balanced_idx].to(device)
    side_effects = [all_se[i] for i in balanced_idx]
    
    print(f"Selected {len(side_effects)} balanced side effects:")
    for i, se in enumerate(side_effects):
        pos_rate = labels[:, i].sum().item() / n_drugs
        print(f"  {se}: {pos_rate:.1%} positive")
    
    in_channels = {nt: data[nt].x.shape[1] for nt in data.node_types}
    
    hgt = HGT(in_channels, 256, 4, 2, data.metadata()).to(device)
    adr_pred = ADRPredictor(256, len(side_effects)).to(device)
    
    optimizer = AdamW(list(hgt.parameters()) + list(adr_pred.parameters()), lr=0.001, weight_decay=0.01)
    
    x_dict = {k: v.to(device) for k, v in data.x_dict.items()}
    edge_index_dict = {k: v.to(device) for k, v in data.edge_index_dict.items()}
    
    perm = torch.randperm(n_drugs)
    train_idx, val_idx = perm[:int(0.8*n_drugs)], perm[int(0.8*n_drugs):]
    
    best_auc = 0
    
    for epoch in range(1, 301):
        hgt.train()
        adr_pred.train()
        
        z = hgt(x_dict, edge_index_dict)
        pred = adr_pred(z['drug'][train_idx])
        loss = F.binary_cross_entropy(pred, labels[train_idx])
        
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        
        if epoch % 10 == 0:
            hgt.eval()
            adr_pred.eval()
            with torch.no_grad():
                z = hgt(x_dict, edge_index_dict)
                pred_val = adr_pred(z['drug'][val_idx]).cpu().numpy()
                labels_np = labels[val_idx].cpu().numpy()
                
                aucs = []
                for i in range(len(side_effects)):
                    if 0 < labels_np[:, i].sum() < len(labels_np):
                        aucs.append(roc_auc_score(labels_np[:, i], pred_val[:, i]))
                avg_auc = np.mean(aucs) if aucs else 0
            
            print(f"Epoch {epoch:03d} | Loss: {loss:.4f} | Val AUC: {avg_auc:.4f}")
            
            if avg_auc > best_auc:
                best_auc = avg_auc
                torch.save({
                    'epoch': epoch,
                    'hgt': hgt.state_dict(),
                    'adr': adr_pred.state_dict(),
                    'auc': avg_auc,
                    'side_effects': side_effects,
                }, DRIVE_CHECKPOINTS / 'best_adr_balanced.pt')
                print(f"  💾 Saved! (AUC: {avg_auc:.4f})")
    
    print(f"\n✅ Best AUC: {best_auc:.4f}")

if __name__ == "__main__":
    train()
