"""Use pretrained HGT (from DTI training) for ADR prediction."""
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

NUM_SE = 20

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
        out = {nt: F.relu(self.lin_dict[nt](x)) for nt, x in x_dict.items()}
        for conv in self.convs:
            out = conv(out, edge_index_dict)
        return out

class ADRClassifier(nn.Module):
    """Larger classifier with residual connections."""
    def __init__(self, in_dim, hidden, num_se):
        super().__init__()
        self.fc1 = nn.Linear(in_dim, hidden)
        self.fc2 = nn.Linear(hidden, hidden)
        self.fc3 = nn.Linear(hidden, hidden)
        self.fc4 = nn.Linear(hidden, num_se)
        self.bn1 = nn.BatchNorm1d(hidden)
        self.bn2 = nn.BatchNorm1d(hidden)
        self.drop = nn.Dropout(0.4)
    
    def forward(self, x):
        x = F.relu(self.bn1(self.fc1(x)))
        res = x
        x = F.relu(self.bn2(self.fc2(x)))
        x = self.drop(x)
        x = F.relu(self.fc3(x)) + res  # Residual
        x = self.drop(x)
        return torch.sigmoid(self.fc4(x))

def train():
    data = torch.load(DRIVE_PROCESSED / 'graph_data.pt', weights_only=False)
    labels_data = torch.load(DRIVE_PROCESSED / 'side_effect_labels.pt', weights_only=False)
    
    labels = labels_data['labels'][:, :NUM_SE].to(device)
    side_effects = labels_data['side_effects'][:NUM_SE]
    
    print(f"Using {NUM_SE} side effects")
    
    in_channels = {nt: data[nt].x.shape[1] for nt in data.node_types}
    
    # Load pretrained HGT from DTI training
    hgt = HGT(in_channels, 256, 4, 2, data.metadata()).to(device)
    
    checkpoint = torch.load(DRIVE_CHECKPOINTS / 'best_model_full.pt', weights_only=False)
    hgt.load_state_dict(checkpoint['model_state_dict'])
    print(f"✅ Loaded pretrained HGT (DTI AUC: {checkpoint['auc_dti']:.4f})")
    
    # Freeze HGT initially
    for param in hgt.parameters():
        param.requires_grad = False
    
    # Classifier with more capacity
    classifier = ADRClassifier(256, 512, NUM_SE).to(device)
    
    x_dict = {k: v.to(device) for k, v in data.x_dict.items()}
    edge_index_dict = {k: v.to(device) for k, v in data.edge_index_dict.items()}
    
    n_drugs = labels.shape[0]
    perm = torch.randperm(n_drugs)
    train_idx, val_idx = perm[:int(0.8*n_drugs)], perm[int(0.8*n_drugs):]
    
    # Phase 1: Train classifier only (50 epochs)
    print("\n📌 Phase 1: Training classifier with frozen HGT...")
    optimizer = AdamW(classifier.parameters(), lr=0.002)
    
    for epoch in range(1, 51):
        classifier.train()
        with torch.no_grad():
            z = hgt(x_dict, edge_index_dict)
        pred = classifier(z['drug'][train_idx])
        loss = F.binary_cross_entropy(pred, labels[train_idx])
        
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        
        if epoch % 10 == 0:
            classifier.eval()
            with torch.no_grad():
                z = hgt(x_dict, edge_index_dict)
                pred_val = classifier(z['drug'][val_idx]).cpu().numpy()
                labels_np = labels[val_idx].cpu().numpy()
                aucs = [roc_auc_score(labels_np[:, i], pred_val[:, i]) 
                        for i in range(NUM_SE) if 0 < labels_np[:, i].sum() < len(labels_np)]
                auc = np.mean(aucs) if aucs else 0
            print(f"Epoch {epoch:03d} | Loss: {loss:.4f} | AUC: {auc:.4f}")
    
    # Phase 2: Fine-tune HGT + classifier together (200 epochs)
    print("\n📌 Phase 2: Fine-tuning HGT + classifier...")
    for param in hgt.parameters():
        param.requires_grad = True
    
    optimizer = AdamW(list(hgt.parameters()) + list(classifier.parameters()), lr=0.0005, weight_decay=0.01)
    
    best_auc = 0
    
    for epoch in range(51, 251):
        hgt.train()
        classifier.train()
        
        z = hgt(x_dict, edge_index_dict)
        pred = classifier(z['drug'][train_idx])
        loss = F.binary_cross_entropy(pred, labels[train_idx])
        
        optimizer.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(list(hgt.parameters()) + list(classifier.parameters()), 1.0)
        optimizer.step()
        
        if epoch % 10 == 0:
            hgt.eval()
            classifier.eval()
            with torch.no_grad():
                z = hgt(x_dict, edge_index_dict)
                pred_val = classifier(z['drug'][val_idx]).cpu().numpy()
                labels_np = labels[val_idx].cpu().numpy()
                aucs = [roc_auc_score(labels_np[:, i], pred_val[:, i]) 
                        for i in range(NUM_SE) if 0 < labels_np[:, i].sum() < len(labels_np)]
                auc = np.mean(aucs) if aucs else 0
            
            print(f"Epoch {epoch:03d} | Loss: {loss:.4f} | AUC: {auc:.4f}")
            
            if auc > best_auc:
                best_auc = auc
                torch.save({
                    'hgt': hgt.state_dict(),
                    'classifier': classifier.state_dict(),
                    'auc': auc,
                    'side_effects': side_effects,
                }, DRIVE_CHECKPOINTS / 'best_adr_final.pt')
                print(f"  💾 Saved! (AUC: {auc:.4f})")
    
    print(f"\n✅ Final Best AUC: {best_auc:.4f}")

if __name__ == "__main__":
    train()
