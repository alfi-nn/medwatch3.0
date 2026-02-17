"""Train ADR (side effect) prediction using HGT embeddings."""
import torch
import torch.nn as nn
import torch.nn.functional as F
from pathlib import Path
from torch.optim import Adam
from sklearn.metrics import roc_auc_score, average_precision_score
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
        out = {nt: self.lin_dict[nt](x).relu() for nt, x in x_dict.items()}
        for conv in self.convs:
            out = conv(out, edge_index_dict)
        return out

class ADRPredictor(nn.Module):
    def __init__(self, hidden_channels, num_side_effects):
        super().__init__()
        self.fc1 = Linear(hidden_channels, hidden_channels)
        self.fc2 = Linear(hidden_channels, num_side_effects)
    
    def forward(self, x):
        x = F.relu(self.fc1(x))
        return torch.sigmoid(self.fc2(x))

def train():
    # Load data
    data = torch.load(DRIVE_PROCESSED / 'graph_data.pt', weights_only=False)
    labels_data = torch.load(DRIVE_PROCESSED / 'side_effect_labels.pt', weights_only=False)
    
    labels = labels_data['labels'].to(device)
    side_effects = labels_data['side_effects']
    num_se = len(side_effects)
    
    print(f"Side effects: {num_se}")
    print(f"Labels shape: {labels.shape}")
    
    in_channels = {nt: data[nt].x.shape[1] for nt in data.node_types}
    
    # Models
    hgt = HGT(in_channels, 256, 4, 2, data.metadata()).to(device)
    adr_pred = ADRPredictor(256, num_se).to(device)
    
    # Load pretrained HGT weights
    checkpoint = torch.load(DRIVE_CHECKPOINTS / 'best_model_full.pt', weights_only=False)
    hgt.load_state_dict(checkpoint['model_state_dict'])
    print("✅ Loaded pretrained HGT weights")
    
    optimizer = Adam(list(hgt.parameters()) + list(adr_pred.parameters()), lr=0.0005)
    
    x_dict = {k: v.to(device) for k, v in data.x_dict.items()}
    edge_index_dict = {k: v.to(device) for k, v in data.edge_index_dict.items()}
    
    # Train/val split
    n_drugs = labels.shape[0]
    perm = torch.randperm(n_drugs)
    train_idx = perm[:int(0.8 * n_drugs)]
    val_idx = perm[int(0.8 * n_drugs):]
    
    best_auc = 0
    
    for epoch in range(1, 101):
        hgt.train()
        adr_pred.train()
        
        z_dict = hgt(x_dict, edge_index_dict)
        drug_emb = z_dict['drug']
        
        pred = adr_pred(drug_emb[train_idx])
        loss = F.binary_cross_entropy(pred, labels[train_idx])
        
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        
        if epoch % 10 == 0:
            hgt.eval()
            adr_pred.eval()
            with torch.no_grad():
                z_dict = hgt(x_dict, edge_index_dict)
                pred_val = adr_pred(z_dict['drug'][val_idx])
                
                pred_np = pred_val.cpu().numpy()
                labels_np = labels[val_idx].cpu().numpy()
                
                # Compute AUC per side effect
                aucs = []
                for i in range(num_se):
                    if labels_np[:, i].sum() > 0:
                        auc = roc_auc_score(labels_np[:, i], pred_np[:, i])
                        aucs.append(auc)
                
                avg_auc = np.mean(aucs) if aucs else 0
            
            print(f"Epoch {epoch:03d} | Loss: {loss:.4f} | Val AUC: {avg_auc:.4f}")
            
            if avg_auc > best_auc:
                best_auc = avg_auc
                torch.save({
                    'epoch': epoch,
                    'hgt_state_dict': hgt.state_dict(),
                    'adr_state_dict': adr_pred.state_dict(),
                    'auc': avg_auc,
                    'side_effects': side_effects,
                }, DRIVE_CHECKPOINTS / 'best_adr_model.pt')
                print(f"  💾 Saved ADR model (AUC: {avg_auc:.4f})")
    
    print(f"\n✅ ADR Training complete! Best AUC: {best_auc:.4f}")

if __name__ == "__main__":
    train()
