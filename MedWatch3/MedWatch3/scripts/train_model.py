"""
Train HGT model with checkpointing to Google Drive.
Fixed dimension handling.
"""
import os
import torch
import torch.nn.functional as F
from pathlib import Path
from torch.optim import Adam
from sklearn.metrics import roc_auc_score
from torch_geometric.nn import HGTConv, Linear

DRIVE_PROCESSED = Path(os.environ.get('DRIVE_DATA_PROCESSED', 'data/processed'))
DRIVE_CHECKPOINTS = Path(os.environ.get('DRIVE_CHECKPOINTS', 'checkpoints'))
device = 'cuda' if torch.cuda.is_available() else 'cpu'

class HGT(torch.nn.Module):
    def __init__(self, in_channels_dict, hidden_channels, num_heads, num_layers, metadata):
        super().__init__()
        self.lin_dict = torch.nn.ModuleDict()
        for node_type in metadata[0]:
            # Explicit input dimension instead of lazy (-1)
            in_dim = in_channels_dict.get(node_type, hidden_channels)
            self.lin_dict[node_type] = Linear(in_dim, hidden_channels)
        
        self.convs = torch.nn.ModuleList()
        for _ in range(num_layers):
            self.convs.append(HGTConv(hidden_channels, hidden_channels, metadata, num_heads))
    
    def forward(self, x_dict, edge_index_dict):
        out_dict = {}
        for node_type, x in x_dict.items():
            out_dict[node_type] = self.lin_dict[node_type](x).relu()
        for conv in self.convs:
            out_dict = conv(out_dict, edge_index_dict)
        return out_dict

class LinkPredictor(torch.nn.Module):
    def __init__(self, hidden_channels):
        super().__init__()
        self.lin1 = Linear(hidden_channels * 2, hidden_channels)
        self.lin2 = Linear(hidden_channels, 1)
    
    def forward(self, x_src, x_dst):
        x = torch.cat([x_src, x_dst], dim=-1)
        return torch.sigmoid(self.lin2(self.lin1(x).relu()))

def train():
    # Load graph
    data = torch.load(DRIVE_PROCESSED / 'graph_data.pt', weights_only=False)
    print(f"Loaded graph: {data}")
    
    # Get input dimensions from data
    in_channels_dict = {node_type: data[node_type].x.shape[1] 
                        for node_type in data.node_types}
    print(f"Input dimensions: {in_channels_dict}")
    
    # Model with explicit input dimensions
    model = HGT(in_channels_dict, 256, 4, 2, data.metadata()).to(device)
    predictor = LinkPredictor(256).to(device)
    optimizer = Adam(list(model.parameters()) + list(predictor.parameters()), lr=0.001)
    
    # Move data to device
    x_dict = {k: v.to(device) for k, v in data.x_dict.items()}
    edge_index_dict = {k: v.to(device) for k, v in data.edge_index_dict.items()}
    edge_index = data['drug', 'interacts', 'drug'].edge_index.to(device)
    
    best_auc = 0
    
    for epoch in range(1, 101):
        model.train()
        predictor.train()
        
        # Forward pass
        z_dict = model(x_dict, edge_index_dict)
        
        # Positive samples
        pos_src, pos_dst = edge_index[0], edge_index[1]
        pos_pred = predictor(z_dict['drug'][pos_src], z_dict['drug'][pos_dst])
        
        # Negative samples
        neg_dst = torch.randint(0, data['drug'].num_nodes, pos_dst.shape, device=device)
        neg_pred = predictor(z_dict['drug'][pos_src], z_dict['drug'][neg_dst])
        
        # Loss
        loss = F.binary_cross_entropy(pos_pred, torch.ones_like(pos_pred)) + \
               F.binary_cross_entropy(neg_pred, torch.zeros_like(neg_pred))
        
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        
        # Evaluate
        if epoch % 10 == 0:
            model.eval()
            with torch.no_grad():
                preds = torch.cat([pos_pred, neg_pred]).cpu().numpy().flatten()
                labels = torch.cat([torch.ones(len(pos_pred)), torch.zeros(len(neg_pred))]).numpy()
                auc = roc_auc_score(labels, preds)
            
            print(f"Epoch {epoch:03d} | Loss: {loss:.4f} | AUC: {auc:.4f}")
            
            # Save checkpoint
            if auc > best_auc:
                best_auc = auc
                torch.save({
                    'epoch': epoch,
                    'model_state_dict': model.state_dict(),
                    'predictor_state_dict': predictor.state_dict(),
                    'auc': auc,
                }, DRIVE_CHECKPOINTS / 'best_model.pt')
                print(f"  💾 Saved checkpoint (AUC: {auc:.4f})")
    
    print(f"\n✅ Training complete! Best AUC: {best_auc:.4f}")

if __name__ == "__main__":
    train()
