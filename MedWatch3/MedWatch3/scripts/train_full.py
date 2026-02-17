"""Train HGT on full graph with drugs + proteins."""
import os
import torch
import torch.nn.functional as F
from pathlib import Path
from torch.optim import Adam
from sklearn.metrics import roc_auc_score
from torch_geometric.nn import HGTConv, Linear

DRIVE_PROCESSED = Path('/content/drive/MyDrive/MedWatch3/data/processed')
DRIVE_CHECKPOINTS = Path('/content/drive/MyDrive/MedWatch3/checkpoints')
device = 'cuda' if torch.cuda.is_available() else 'cpu'

class HGT(torch.nn.Module):
    def __init__(self, in_channels_dict, hidden_channels, num_heads, num_layers, metadata):
        super().__init__()
        self.lin_dict = torch.nn.ModuleDict()
        for node_type in metadata[0]:
            in_dim = in_channels_dict[node_type]
            self.lin_dict[node_type] = Linear(in_dim, hidden_channels)
        
        self.convs = torch.nn.ModuleList()
        for _ in range(num_layers):
            self.convs.append(HGTConv(hidden_channels, hidden_channels, metadata, num_heads))
    
    def forward(self, x_dict, edge_index_dict):
        out = {}
        for node_type, x in x_dict.items():
            out[node_type] = self.lin_dict[node_type](x).relu()
        for conv in self.convs:
            out = conv(out, edge_index_dict)
        return out

class LinkPredictor(torch.nn.Module):
    def __init__(self, hidden_channels):
        super().__init__()
        self.lin1 = Linear(hidden_channels * 2, hidden_channels)
        self.lin2 = Linear(hidden_channels, 1)
    
    def forward(self, x_src, x_dst):
        x = torch.cat([x_src, x_dst], dim=-1)
        return torch.sigmoid(self.lin2(self.lin1(x).relu()))

def train():
    data = torch.load(DRIVE_PROCESSED / 'graph_data.pt', weights_only=False)
    print(f"Graph: {data}")
    
    in_channels = {nt: data[nt].x.shape[1] for nt in data.node_types}
    print(f"Input dims: {in_channels}")
    
    model = HGT(in_channels, 256, 4, 2, data.metadata()).to(device)
    predictor_ddi = LinkPredictor(256).to(device)
    predictor_dti = LinkPredictor(256).to(device)
    
    optimizer = Adam(
        list(model.parameters()) + list(predictor_ddi.parameters()) + list(predictor_dti.parameters()),
        lr=0.001
    )
    
    x_dict = {k: v.to(device) for k, v in data.x_dict.items()}
    edge_index_dict = {k: v.to(device) for k, v in data.edge_index_dict.items()}
    
    ddi_edge = data['drug', 'interacts', 'drug'].edge_index.to(device)
    dti_edge = data['drug', 'targets', 'protein'].edge_index.to(device)
    
    best_auc = 0
    
    for epoch in range(1, 151):
        model.train()
        
        z_dict = model(x_dict, edge_index_dict)
        
        # DDI loss
        pos_ddi = predictor_ddi(z_dict['drug'][ddi_edge[0]], z_dict['drug'][ddi_edge[1]])
        neg_dst = torch.randint(0, data['drug'].num_nodes, ddi_edge[1].shape, device=device)
        neg_ddi = predictor_ddi(z_dict['drug'][ddi_edge[0]], z_dict['drug'][neg_dst])
        loss_ddi = F.binary_cross_entropy(pos_ddi, torch.ones_like(pos_ddi)) + \
                   F.binary_cross_entropy(neg_ddi, torch.zeros_like(neg_ddi))
        
        # DTI loss
        pos_dti = predictor_dti(z_dict['drug'][dti_edge[0]], z_dict['protein'][dti_edge[1]])
        neg_prot = torch.randint(0, data['protein'].num_nodes, dti_edge[1].shape, device=device)
        neg_dti = predictor_dti(z_dict['drug'][dti_edge[0]], z_dict['protein'][neg_prot])
        loss_dti = F.binary_cross_entropy(pos_dti, torch.ones_like(pos_dti)) + \
                   F.binary_cross_entropy(neg_dti, torch.zeros_like(neg_dti))
        
        loss = loss_ddi + loss_dti
        
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        
        if epoch % 10 == 0:
            model.eval()
            with torch.no_grad():
                # DDI AUC
                preds = torch.cat([pos_ddi, neg_ddi]).cpu().numpy().flatten()
                labels = torch.cat([torch.ones(len(pos_ddi)), torch.zeros(len(neg_ddi))]).numpy()
                auc_ddi = roc_auc_score(labels, preds)
                
                # DTI AUC
                preds = torch.cat([pos_dti, neg_dti]).cpu().numpy().flatten()
                labels = torch.cat([torch.ones(len(pos_dti)), torch.zeros(len(neg_dti))]).numpy()
                auc_dti = roc_auc_score(labels, preds)
            
            print(f"Epoch {epoch:03d} | Loss: {loss:.4f} | DDI AUC: {auc_ddi:.4f} | DTI AUC: {auc_dti:.4f}")
            
            avg_auc = (auc_ddi + auc_dti) / 2
            if avg_auc > best_auc:
                best_auc = avg_auc
                torch.save({
                    'epoch': epoch,
                    'model_state_dict': model.state_dict(),
                    'auc_ddi': auc_ddi,
                    'auc_dti': auc_dti,
                }, DRIVE_CHECKPOINTS / 'best_model_full.pt')
                print(f"  💾 Saved! (Avg AUC: {avg_auc:.4f})")
    
    print(f"\n✅ Training complete! Best Avg AUC: {best_auc:.4f}")

if __name__ == "__main__":
    train()
