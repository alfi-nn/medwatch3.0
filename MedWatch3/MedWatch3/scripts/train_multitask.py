"""Multi-task learning: Joint DDI + DTI + ADR training."""
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

NUM_SE = 20
EPOCHS = 300
LR = 0.001

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
            out = {k: F.dropout(v, 0.1, self.training) for k, v in out.items()}
        return out

class LinkPredictor(nn.Module):
    def __init__(self, hidden_channels):
        super().__init__()
        self.lin1 = Linear(hidden_channels * 2, hidden_channels)
        self.lin2 = Linear(hidden_channels, 1)
    
    def forward(self, x_src, x_dst):
        x = torch.cat([x_src, x_dst], dim=-1)
        return torch.sigmoid(self.lin2(F.relu(self.lin1(x))))

class ADRPredictor(nn.Module):
    def __init__(self, hidden_channels, num_se):
        super().__init__()
        self.fc1 = Linear(hidden_channels, hidden_channels * 2)
        self.fc2 = Linear(hidden_channels * 2, hidden_channels)
        self.fc3 = Linear(hidden_channels, num_se)
        self.drop = nn.Dropout(0.3)
    
    def forward(self, x):
        x = self.drop(F.relu(self.fc1(x)))
        x = self.drop(F.relu(self.fc2(x)))
        return torch.sigmoid(self.fc3(x))

def train():
    # Load data
    data = torch.load(DRIVE_PROCESSED / 'graph_data.pt', weights_only=False)
    labels_data = torch.load(DRIVE_PROCESSED / 'side_effect_labels.pt', weights_only=False)
    
    labels = labels_data['labels'][:, :NUM_SE].to(device)
    side_effects = labels_data['side_effects'][:NUM_SE]
    
    print(f"Graph: {data}")
    print(f"Side effects: {NUM_SE}")
    
    in_channels = {nt: data[nt].x.shape[1] for nt in data.node_types}
    
    # Models
    hgt = HGT(in_channels, 256, 4, 2, data.metadata()).to(device)
    ddi_pred = LinkPredictor(256).to(device)
    dti_pred = LinkPredictor(256).to(device)
    adr_pred = ADRPredictor(256, NUM_SE).to(device)
    
    all_params = (list(hgt.parameters()) + list(ddi_pred.parameters()) + 
                  list(dti_pred.parameters()) + list(adr_pred.parameters()))
    
    optimizer = AdamW(all_params, lr=LR, weight_decay=0.01)
    scheduler = CosineAnnealingLR(optimizer, T_max=EPOCHS)
    
    x_dict = {k: v.to(device) for k, v in data.x_dict.items()}
    edge_index_dict = {k: v.to(device) for k, v in data.edge_index_dict.items()}
    
    ddi_edge = data['drug', 'interacts', 'drug'].edge_index.to(device)
    dti_edge = data['drug', 'targets', 'protein'].edge_index.to(device)
    
    # ADR train/val split
    n_drugs = labels.shape[0]
    perm = torch.randperm(n_drugs)
    train_idx = perm[:int(0.8 * n_drugs)]
    val_idx = perm[int(0.8 * n_drugs):]
    
    best_combined = 0
    
    for epoch in range(1, EPOCHS + 1):
        hgt.train()
        ddi_pred.train()
        dti_pred.train()
        adr_pred.train()
        
        z_dict = hgt(x_dict, edge_index_dict)
        
        # === DDI Loss ===
        pos_ddi = ddi_pred(z_dict['drug'][ddi_edge[0]], z_dict['drug'][ddi_edge[1]])
        neg_dst = torch.randint(0, data['drug'].num_nodes, ddi_edge[1].shape, device=device)
        neg_ddi = ddi_pred(z_dict['drug'][ddi_edge[0]], z_dict['drug'][neg_dst])
        loss_ddi = F.binary_cross_entropy(pos_ddi, torch.ones_like(pos_ddi)) + \
                   F.binary_cross_entropy(neg_ddi, torch.zeros_like(neg_ddi))
        
        # === DTI Loss ===
        pos_dti = dti_pred(z_dict['drug'][dti_edge[0]], z_dict['protein'][dti_edge[1]])
        neg_prot = torch.randint(0, data['protein'].num_nodes, dti_edge[1].shape, device=device)
        neg_dti = dti_pred(z_dict['drug'][dti_edge[0]], z_dict['protein'][neg_prot])
        loss_dti = F.binary_cross_entropy(pos_dti, torch.ones_like(pos_dti)) + \
                   F.binary_cross_entropy(neg_dti, torch.zeros_like(neg_dti))
        
        # === ADR Loss (weighted higher) ===
        adr_out = adr_pred(z_dict['drug'][train_idx])
        loss_adr = F.binary_cross_entropy(adr_out, labels[train_idx])
        
        # Combined loss (weight ADR higher)
        loss = loss_ddi + loss_dti + 3.0 * loss_adr
        
        optimizer.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(all_params, 1.0)
        optimizer.step()
        scheduler.step()
        
        if epoch % 10 == 0:
            hgt.eval()
            with torch.no_grad():
                z_dict = hgt(x_dict, edge_index_dict)
                
                # DDI AUC
                preds = torch.cat([pos_ddi, neg_ddi]).cpu().numpy().flatten()
                labs = np.concatenate([np.ones(len(pos_ddi)), np.zeros(len(neg_ddi))])
                auc_ddi = roc_auc_score(labs, preds)
                
                # DTI AUC
                preds = torch.cat([pos_dti, neg_dti]).cpu().numpy().flatten()
                labs = np.concatenate([np.ones(len(pos_dti)), np.zeros(len(neg_dti))])
                auc_dti = roc_auc_score(labs, preds)
                
                # ADR AUC
                pred_val = adr_pred(z_dict['drug'][val_idx]).cpu().numpy()
                labels_np = labels[val_idx].cpu().numpy()
                aucs = [roc_auc_score(labels_np[:, i], pred_val[:, i]) 
                        for i in range(NUM_SE) 
                        if 0 < labels_np[:, i].sum() < len(labels_np)]
                auc_adr = np.mean(aucs) if aucs else 0
            
            combined = (auc_ddi + auc_dti + auc_adr) / 3
            print(f"Epoch {epoch:03d} | Loss: {loss:.4f} | DDI: {auc_ddi:.4f} | DTI: {auc_dti:.4f} | ADR: {auc_adr:.4f}")
            
            if combined > best_combined:
                best_combined = combined
                torch.save({
                    'epoch': epoch,
                    'hgt': hgt.state_dict(),
                    'ddi_pred': ddi_pred.state_dict(),
                    'dti_pred': dti_pred.state_dict(),
                    'adr_pred': adr_pred.state_dict(),
                    'auc_ddi': auc_ddi,
                    'auc_dti': auc_dti,
                    'auc_adr': auc_adr,
                    'side_effects': side_effects,
                }, DRIVE_CHECKPOINTS / 'best_multitask.pt')
                print(f"  💾 Saved! (Combined: {combined:.4f})")
    
    print(f"\n✅ Multi-task training complete!")

if __name__ == "__main__":
    train()
