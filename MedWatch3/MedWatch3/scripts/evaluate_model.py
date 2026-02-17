"""
Comprehensive model evaluation: AUC, F1, precision, recall, calibration, 
per-class breakdown, and drug-discrimination analysis.
"""
import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from pathlib import Path
from sklearn.metrics import (
    roc_auc_score, f1_score, precision_score, recall_score,
    accuracy_score, average_precision_score
)
from torch_geometric.nn import HGTConv, Linear
import sys

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / 'data' / 'processed'
CHECKPOINT_DIR = BASE_DIR / 'checkpoints'
OUTPUT = BASE_DIR / 'evaluation_report.txt'
DEVICE = 'cpu'

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

class ADRPredictor3Layer(nn.Module):
    def __init__(self, hidden, num_se):
        super().__init__()
        self.fc1 = Linear(hidden, hidden * 2)
        self.fc2 = Linear(hidden * 2, hidden)
        self.fc3 = Linear(hidden, num_se)
        self.dropout = nn.Dropout(0.3)
    
    def forward(self, x):
        x = F.relu(self.fc1(x))
        x = self.dropout(x)
        x = F.relu(self.fc2(x))
        x = self.dropout(x)
        return self.fc3(x)

def log(msg, f):
    print(msg)
    f.write(msg + '\n')

def evaluate_checkpoint(ckpt_name, data, labels, side_effects, f):
    """Evaluate a single checkpoint."""
    ckpt_path = CHECKPOINT_DIR / ckpt_name
    if not ckpt_path.exists():
        log(f"  Checkpoint {ckpt_name} not found, skipping.", f)
        return None
    
    ckpt = torch.load(ckpt_path, weights_only=False, map_location='cpu')
    
    # Determine number of side effects in checkpoint
    if 'side_effects' in ckpt:
        ckpt_se = ckpt['side_effects']
        num_se = len(ckpt_se)
    else:
        # Infer from state dict
        for k, v in ckpt.get('adr_state_dict', {}).items():
            if 'fc3' in k and 'weight' in k:
                num_se = v.shape[0]
                break
        ckpt_se = side_effects[:num_se] if num_se <= len(side_effects) else side_effects
    
    log(f"\n  Checkpoint: {ckpt_name}", f)
    log(f"  Side effects: {num_se}", f)
    if 'auc' in ckpt:
        log(f"  Saved AUC: {ckpt['auc']:.4f}", f)
    if 'epoch' in ckpt:
        log(f"  Trained epochs: {ckpt['epoch']}", f)
    
    # Load HGT
    in_channels = {nt: data[nt].x.shape[1] for nt in data.node_types}
    hgt = HGT(in_channels, 256, 4, 2, data.metadata()).to(DEVICE)
    
    if 'hgt_state_dict' in ckpt:
        hgt.load_state_dict(ckpt['hgt_state_dict'])
    elif 'model_state_dict' in ckpt:
        try:
            hgt.load_state_dict(ckpt['model_state_dict'])
        except:
            log(f"  Warning: Could not load HGT state dict", f)
    
    hgt.eval()
    
    # Load ADR predictor
    adr = ADRPredictor3Layer(256, num_se).to(DEVICE)
    if 'adr_state_dict' in ckpt:
        try:
            adr.load_state_dict(ckpt['adr_state_dict'])
        except Exception as e:
            log(f"  Warning: ADR load error: {e}", f)
            return None
    adr.eval()
    
    # Get predictions on original drugs
    x_dict = {k: v.to(DEVICE) for k, v in data.x_dict.items()}
    edge_index_dict = {k: v.to(DEVICE) for k, v in data.edge_index_dict.items()}
    
    with torch.no_grad():
        z_dict = hgt(x_dict, edge_index_dict)
        drug_emb = z_dict['drug'][:labels.shape[0]]  # Only original drugs
        logits = adr(drug_emb)
        probs = torch.sigmoid(logits).cpu().numpy()
    
    # Use matching labels
    use_labels = labels[:, :num_se].numpy() if labels.shape[1] >= num_se else labels.numpy()
    if use_labels.shape[1] != probs.shape[1]:
        log(f"  Skipping: label/pred shape mismatch ({use_labels.shape} vs {probs.shape})", f)
        return None
    
    # Global metrics
    preds_binary = (probs > 0.5).astype(float)
    
    # Per-class AUC
    aucs = []
    aps = []
    f1s = []
    for i in range(num_se):
        col = use_labels[:, i]
        if 0 < col.sum() < len(col):
            try:
                auc = roc_auc_score(col, probs[:, i])
                aucs.append(auc)
            except:
                pass
            try:
                ap = average_precision_score(col, probs[:, i])
                aps.append(ap)
            except:
                pass
        f1 = f1_score(col, preds_binary[:, i], zero_division=0)
        f1s.append(f1)
    
    avg_auc = np.mean(aucs) if aucs else 0
    avg_ap = np.mean(aps) if aps else 0
    avg_f1 = np.mean(f1s)
    
    # Calibration: how close are average predictions to average labels?
    mean_pred = probs.mean()
    mean_label = use_labels.mean()
    calibration_error = abs(mean_pred - mean_label)
    
    # Discrimination: std of predictions across drugs (higher = better separation)
    pred_std_per_se = probs.std(axis=0).mean()
    pred_range = probs.max() - probs.min()
    
    log(f"  Avg AUC: {avg_auc:.4f}", f)
    log(f"  Avg AP: {avg_ap:.4f}", f)
    log(f"  Avg F1: {avg_f1:.4f}", f)
    log(f"  Mean prediction: {mean_pred:.4f} (label mean: {mean_label:.4f})", f)
    log(f"  Calibration error: {calibration_error:.4f}", f)
    log(f"  Pred std across drugs (discrimination): {pred_std_per_se:.4f}", f)
    log(f"  Pred range: {probs.min():.4f} - {probs.max():.4f}", f)
    
    # Top 5 best and worst AUC side effects
    if aucs and len(ckpt_se) == len(aucs):
        auc_pairs = list(zip(ckpt_se, aucs))
        auc_pairs.sort(key=lambda x: x[1], reverse=True)
        log(f"\n  Top 5 best-predicted side effects:", f)
        for name, auc in auc_pairs[:5]:
            log(f"    {name}: AUC={auc:.4f}", f)
        log(f"  Top 5 worst-predicted side effects:", f)
        for name, auc in auc_pairs[-5:]:
            log(f"    {name}: AUC={auc:.4f}", f)
    
    return {
        'auc': avg_auc, 'ap': avg_ap, 'f1': avg_f1,
        'calibration': calibration_error, 'discrimination': pred_std_per_se,
        'range': (probs.min(), probs.max()),
    }

def main():
    with open(OUTPUT, 'w', encoding='utf-8') as f:
        log("=" * 70, f)
        log("COMPREHENSIVE MODEL EVALUATION REPORT", f)
        log("=" * 70, f)
        
        # Load original data 
        data_orig = torch.load(DATA_DIR / 'graph_data.pt', weights_only=False, map_location='cpu')
        labels_orig = torch.load(DATA_DIR / 'side_effect_labels.pt', weights_only=False, map_location='cpu')
        labels_matrix = labels_orig['labels']
        side_effects = labels_orig['side_effects']
        
        log(f"\nDataset: {labels_matrix.shape[0]} drugs, {len(side_effects)} side effects", f)
        
        # Load augmented data for balanced model
        data_aug = torch.load(DATA_DIR / 'graph_data_augmented.pt', weights_only=False, map_location='cpu')
        labels_aug = torch.load(DATA_DIR / 'side_effect_labels_augmented.pt', weights_only=False, map_location='cpu')
        labels_aug_matrix = labels_aug['labels']
        side_effects_aug = labels_aug['side_effects']
        
        # Evaluate each checkpoint
        checkpoints = [
            ('best_adr_v3.pt', data_aug, labels_aug_matrix, side_effects_aug),
            ('best_adr_balanced_v2.pt', data_aug, labels_aug_matrix, side_effects_aug),
            ('best_adr_improved(0.65).pt', data_orig, labels_matrix, side_effects),
            ('best_adr_model.pt', data_orig, labels_matrix, side_effects),
        ]
        
        results = {}
        for ckpt_name, data, labels, se in checkpoints:
            log(f"\n{'─' * 50}", f)
            r = evaluate_checkpoint(ckpt_name, data, labels, se, f)
            if r:
                results[ckpt_name] = r
        
        # Summary comparison table
        log(f"\n{'=' * 70}", f)
        log("SUMMARY COMPARISON", f)
        log(f"{'=' * 70}", f)
        log(f"{'Model':<35} {'AUC':>7} {'AP':>7} {'F1':>7} {'CalErr':>7} {'Discr':>7} {'Range':>15}", f)
        log("-" * 90, f)
        for name, r in results.items():
            rng = f"{r['range'][0]:.2f}-{r['range'][1]:.2f}"
            log(f"{name:<35} {r['auc']:>7.4f} {r['ap']:>7.4f} {r['f1']:>7.4f} {r['calibration']:>7.4f} {r['discrimination']:>7.4f} {rng:>15}", f)
        
        log(f"\n{'=' * 70}", f)
        log("ANALYSIS & RECOMMENDATIONS", f)
        log(f"{'=' * 70}", f)
        
        # Analyze issues
        if results:
            best = max(results.items(), key=lambda x: x[1]['auc'])
            log(f"\nBest model by AUC: {best[0]} (AUC={best[1]['auc']:.4f})", f)
            
            log("\nKey Issues:", f)
            for name, r in results.items():
                if r['auc'] < 0.65:
                    log(f"  [ISSUE] {name}: AUC {r['auc']:.4f} is below 0.65 (random=0.50)", f)
                if r['discrimination'] < 0.05:
                    log(f"  [ISSUE] {name}: Low discrimination ({r['discrimination']:.4f}) — predictions are too similar across drugs", f)
                if r['calibration'] > 0.10:
                    log(f"  [ISSUE] {name}: Poor calibration (error={r['calibration']:.4f})", f)
            
            log("\nRecommendations:", f)
            log("  1. INCREASE TRAINING DATA: 500 drugs is very small. Consider using DrugBank or ChEMBL for more drugs.", f)
            log("  2. USE PRETRAINED HGT: Load HGT weights from DTI task (best_model_full.pt) instead of training from scratch.", f)
            log("  3. ADD NEGATIVE SAMPLING: Instead of random augmentation, use structure-aware negatives.", f)
            log("  4. MULTI-TASK LEARNING: Train side effect prediction jointly with DTI to share representations.", f)
            log("  5. CROSS-VALIDATION: Use 5-fold CV instead of single 80/20 split for more reliable metrics.", f)

if __name__ == "__main__":
    main()
