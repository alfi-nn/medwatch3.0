"""
Updated prediction script with Ensemble (v3 + balanced_v2) and Temperature Scaling.
"""
import torch
import torch.nn as nn
import torch.nn.functional as F
from pathlib import Path
from transformers import AutoTokenizer, AutoModel
import sys

BASE_DIR = Path(__file__).resolve().parent.parent
sys.path.append(str(BASE_DIR))

# Constants
CHECKPOINT_DIR = BASE_DIR / 'checkpoints'
DATA_DIR = BASE_DIR / 'data' / 'processed'
device = 'cuda' if torch.cuda.is_available() else 'cpu'

# Config
TEMPERATURE = 3.5  # Soften predictions (higher = softer, spreads probabilities more)

# Global models
tokenizer = None
chemberta = None
hgt_model = None
adr_predictors = []  # List of (model, weight) tuples for ensemble

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
        out = {nt: F.dropout(self.lin_dict[nt](x).relu(), p=0.1, training=False) 
               for nt, x in x_dict.items()}
        for conv in self.convs:
            out = conv(out, edge_index_dict)
        return out

from torch_geometric.nn import HGTConv, Linear

class ADRPredictor(nn.Module):
    def __init__(self, hidden, num_se, style='v3'):
        super().__init__()
        self.style = style
        if style == 'v3':
            # V3 Architecture with BN
            self.fc1 = Linear(hidden, hidden * 2)
            self.bn1 = nn.BatchNorm1d(hidden * 2)
            self.fc2 = Linear(hidden * 2, hidden)
            self.bn2 = nn.BatchNorm1d(hidden)
            self.fc3 = Linear(hidden, num_se)
            self.dropout = nn.Dropout(0.4)
        else:
            # V2 / Improved Architecture (Sequential)
            self.net = nn.Sequential(
                Linear(hidden, hidden * 2), nn.ReLU(), nn.Dropout(0.3),
                Linear(hidden * 2, hidden), nn.ReLU(), nn.Dropout(0.3),
                Linear(hidden, num_se) # Raw logits, sigmoid applied later
            )

    def forward(self, x):
        if self.style == 'v3':
            x = F.relu(self.bn1(self.fc1(x)))
            x = self.dropout(x)
            x = F.relu(self.bn2(self.fc2(x)))
            x = self.dropout(x)
            return self.fc3(x) # Return LOGITS
        else:
            return self.net(x) # Return LOGITS

def load_models():
    """Load ChemBERTa, HGT, and Ensemble of ADR models."""
    global tokenizer, chemberta, hgt_model, adr_predictors, side_effects
    
    if tokenizer is not None:
        return

    print("Loading models...")
    
    # 1. ChemBERTa
    print("   Loading ChemBERTa...")
    # Use 77M-MLM which has 384 hidden size (matches HGT input)
    tokenizer = AutoTokenizer.from_pretrained("DeepChem/ChemBERTa-77M-MLM")
    chemberta = AutoModel.from_pretrained("DeepChem/ChemBERTa-77M-MLM").to(device)
    chemberta.eval()
    
    # 2. HGT (Graph backbone)
    print("   Loading HGT model...")
    graph_data = torch.load(DATA_DIR / 'graph_data.pt', weights_only=False, map_location='cpu')
    metadata = graph_data.metadata()
    in_channels = {nt: graph_data[nt].x.shape[1] for nt in graph_data.node_types}
    
    hgt_model = HGT(in_channels, 256, 4, 2, metadata).to(device)
    # Load HGT weights from DTI (best general representation)
    dti_ckpt_path = CHECKPOINT_DIR / 'best_model_full.pt'
    if dti_ckpt_path.exists():
        dti_ckpt = torch.load(dti_ckpt_path, weights_only=False, map_location='cpu')
        hgt_model.load_state_dict(dti_ckpt['model_state_dict'])
    hgt_model.eval()
    
    # 3. ADR Ensemble
    print("   Loading ADR Ensemble...")
    
    # Model 1: V3 (High performance, sharp predictions)
    v3_path = CHECKPOINT_DIR / 'best_adr_v3.pt'
    if v3_path.exists():
        print(f"      Loading {v3_path.name}...")
        ckpt = torch.load(v3_path, weights_only=False, map_location=device)
        se_v3 = ckpt['side_effects']
        model_v3 = ADRPredictor(256, len(se_v3), style='v3').to(device)
        model_v3.load_state_dict(ckpt['adr_state_dict'])
        model_v3.eval()
        adr_predictors.append({'model': model_v3, 'se': se_v3, 'weight': 0.6}) # Higher weight for V3
    
    # Model 2: Balanced V2 (Better calibration, wider range)
    v2_path = CHECKPOINT_DIR / 'best_adr_balanced_v2.pt'
    if v2_path.exists():
        print(f"      Loading {v2_path.name}...")
        ckpt = torch.load(v2_path, weights_only=False, map_location=device)
        se_v2 = ckpt['side_effects']
        
        # Determine num_layers architecture based on file
        # Balanced v2 used 3 layers but Sequential style
        model_v2 = ADRPredictor(256, len(se_v2), style='v2').to(device)
        
        # Load state dict with mapping if needed
        sd = ckpt['adr_state_dict']
        if 'fc1.weight' in sd: # Map to Sequential
             from collections import OrderedDict
             new_sd = OrderedDict()
             for k,v in sd.items():
                 if 'fc1' in k: new_sd[k.replace('fc1', 'net.0')] = v
                 elif 'fc2' in k: new_sd[k.replace('fc2', 'net.3')] = v
                 elif 'fc3' in k: new_sd[k.replace('fc3', 'net.6')] = v
                 else: new_sd[k] = v
             model_v2.load_state_dict(new_sd)
        else:
             model_v2.load_state_dict(sd)
             
        model_v2.eval()
        adr_predictors.append({'model': model_v2, 'se': se_v2, 'weight': 0.4})
    
    if not adr_predictors:
        raise FileNotFoundError("No ADR checkpoints found!")
    
    print(f"   Ensemble loaded with {len(adr_predictors)} models.")

def smiles_to_embedding(smiles):
    inputs = tokenizer(smiles, return_tensors='pt', padding=True, truncation=True, max_length=512).to(device)
    with torch.no_grad():
        outputs = chemberta(**inputs)
        embedding = outputs.last_hidden_state[:, 0, :]
    return embedding

def predict_side_effects(smiles):
    load_models()
    
    drug_emb = smiles_to_embedding(smiles)
    
    # Project through HGT drug linear layer
    # Note: In hetero graphs, HGT projects input features to hidden space first
    drug_proj = hgt_model.lin_dict['drug'](drug_emb)
    
    # Get predictions from each ensemble member
    all_results = {} # key: side_effect, value: weighted_prob_sum
    total_weight = 0
    
    with torch.no_grad():
        for item in adr_predictors:
            model = item['model']
            side_effects = item['se']
            weight = item['weight']
            
            # Get logits
            logits = model(drug_proj)
            
            # Apply Temperature Scaling
            # T > 1 softens the distribution (pushes probs towards 0.5)
            scaled_logits = logits / TEMPERATURE
            probs = torch.sigmoid(scaled_logits).cpu().numpy()[0]
            
            # Aggregate by side effect name
            for se_name, prob in zip(side_effects, probs):
                if se_name not in all_results:
                    all_results[se_name] = 0.0
                all_results[se_name] += prob * weight
            
            # Note: We normalize by sum of weights for each SE individually later 
            # (since different models might have different SE sets)
            # But here our V3 and V2 have SAME filtered SE set (76).
            # If they differed, we'd need to track weight counts per SE.
            # Assuming they match for now (both trained on balanced set).
            
    # Normalize
    # If sets are identical, total_weight is just sum of model weights
    total_model_weight = sum(item['weight'] for item in adr_predictors)
    
    final_preds = {}
    for se, weighted_sum in all_results.items():
        # Check which models have this SE to normalize correctly
        # Optimization: assume intersection is 100% since both used filtered set
        # Cast to standard float for JSON serialization
        final_preds[se] = float(weighted_sum / total_model_weight)
        
    return final_preds

def print_predictions(smiles, results):
    print(f"\nInput SMILES: {smiles}")
    print(f"Ensemble Temperature: {TEMPERATURE}")
    print("\nPREDICTED SIDE EFFECTS")
    print("──────────────────────────────────────────────────")
    
    sorted_res = sorted(results.items(), key=lambda x: x[1], reverse=True)
    
    categories = {'HIGH': [], 'MED': [], 'LOW': []}
    
    for se, prob in sorted_res:
        if prob > 0.90: categories['HIGH'].append((se, prob))
        elif prob > 0.70: categories['MED'].append((se, prob))
        else: categories['LOW'].append((se, prob))
    
    print("\nHIGH RISK (>90%):")
    for se, p in categories['HIGH']: print(f"   {se:<25} {p:.1%} {'█'*20}")
    
    print("\nMEDIUM RISK (70-90%):")
    for se, p in categories['MED']: print(f"   {se:<25} {p:.1%} {'█'*int(p*20)}")
    
    print("\nLOWER RISK (<70%):")
    for se, p in categories['LOW']: 
        bar = '█'*int(p*20)
        print(f"   {se:<25} {p:.1%} {bar}")

if __name__ == "__main__":
    while True:
        s = input("\nSMILES: ").strip()
        if s.lower() == 'quit': break
        try:
            res = predict_side_effects(s)
            print_predictions(s, res)
        except Exception as e:
            print(f"Error: {e}")
