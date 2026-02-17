"""Test trained models - predict side effects for drugs."""
import torch
import torch.nn.functional as F
import pandas as pd
from pathlib import Path
from torch_geometric.nn import HGTConv, Linear
import torch.nn as nn

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
        out = {nt: F.relu(self.lin_dict[nt](x)) for nt, x in x_dict.items()}
        for conv in self.convs:
            out = conv(out, edge_index_dict)
        return out

class ADRPredictor(nn.Module):
    def __init__(self, hidden, num_se):
        super().__init__()
        self.fc1 = Linear(hidden, hidden * 2)
        self.fc2 = Linear(hidden * 2, hidden)
        self.fc3 = Linear(hidden, num_se)
        self.drop = nn.Dropout(0.3)
    
    def forward(self, x):
        x = self.drop(F.relu(self.fc1(x)))
        x = self.drop(F.relu(self.fc2(x)))
        return torch.sigmoid(self.fc3(x))

def load_models():
    """Load trained models."""
    print("Loading models...")
    
    # Load graph
    data = torch.load(DRIVE_PROCESSED / 'graph_data.pt', weights_only=False)
    in_channels = {nt: data[nt].x.shape[1] for nt in data.node_types}
    
    # Load HGT
    hgt = HGT(in_channels, 256, 4, 2, data.metadata()).to(device)
    dti_ckpt = torch.load(DRIVE_CHECKPOINTS / 'best_model_full.pt', weights_only=False)
    hgt.load_state_dict(dti_ckpt['model_state_dict'])
    hgt.eval()
    print(f"✅ Loaded HGT (DTI AUC: {dti_ckpt['auc_dti']:.4f})")
    
    # Load ADR predictor
    adr_ckpt = torch.load(DRIVE_CHECKPOINTS / 'best_adr_improved.pt', weights_only=False)
    num_se = len(adr_ckpt['side_effects'])
    adr_pred = ADRPredictor(256, num_se).to(device)
    adr_pred.load_state_dict(adr_ckpt['adr_state_dict'])
    adr_pred.eval()
    print(f"✅ Loaded ADR predictor (AUC: {adr_ckpt['auc']:.4f})")
    
    return hgt, adr_pred, data, adr_ckpt['side_effects']

def predict_for_drug(drug_idx, hgt, adr_pred, data, side_effects, drug_names):
    """Predict side effects for a specific drug."""
    x_dict = {k: v.to(device) for k, v in data.x_dict.items()}
    edge_index_dict = {k: v.to(device) for k, v in data.edge_index_dict.items()}
    
    with torch.no_grad():
        z_dict = hgt(x_dict, edge_index_dict)
        drug_emb = z_dict['drug'][drug_idx:drug_idx+1]
        predictions = adr_pred(drug_emb).cpu().numpy()[0]
    
    return predictions

def main():
    # Load models
    hgt, adr_pred, data, side_effects = load_models()
    
    # Load drug names
    df_drugs = pd.read_csv(DRIVE_PROCESSED / 'nodes_drugs.csv')
    drug_names = df_drugs['drug_name'].tolist() if 'drug_name' in df_drugs.columns else df_drugs['drug_id'].tolist()
    
    print(f"\n📊 Testing on {len(drug_names)} drugs...")
    print("=" * 60)
    
    # Test on random drugs
    import random
    test_indices = random.sample(range(len(drug_names)), min(5, len(drug_names)))
    
    for idx in test_indices:
        drug_name = drug_names[idx]
        predictions = predict_for_drug(idx, hgt, adr_pred, data, side_effects, drug_names)
        
        print(f"\n💊 Drug: {drug_name}")
        print("-" * 40)
        
        # Sort by probability
        sorted_idx = predictions.argsort()[::-1]
        
        print("Top 5 predicted side effects:")
        for i, se_idx in enumerate(sorted_idx[:5]):
            prob = predictions[se_idx]
            se_name = side_effects[se_idx]
            bar = "█" * int(prob * 20)
            print(f"  {i+1}. {se_name}: {prob:.1%} {bar}")
        
        print("\nLeast likely side effects:")
        for i, se_idx in enumerate(sorted_idx[-3:]):
            prob = predictions[se_idx]
            se_name = side_effects[se_idx]
            print(f"  • {se_name}: {prob:.1%}")
    
    print("\n" + "=" * 60)
    print("✅ Testing complete!")

if __name__ == "__main__":
    main()
