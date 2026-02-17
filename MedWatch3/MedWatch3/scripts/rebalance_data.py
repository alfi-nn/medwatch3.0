"""
Phase 1: Rebalance side effect labels by filtering to discriminative side effects.
Phase 2: Augment graph with synthetic drug nodes for better negative coverage.
"""
import torch
import numpy as np
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / 'data' / 'processed'
OUTPUT_FILE = BASE_DIR / 'rebalance_log.txt'

def log(msg, f=None):
    print(msg)
    if f:
        f.write(msg + '\n')

def rebalance_labels():
    """Filter side effects to keep only those with 20-80% positive rate."""
    with open(OUTPUT_FILE, 'w', encoding='utf-8') as f:
        log("=" * 60, f)
        log("PHASE 1: REBALANCING SIDE EFFECT LABELS", f)
        log("=" * 60, f)
        
        # Load original labels
        labels_data = torch.load(DATA_DIR / 'side_effect_labels.pt', weights_only=False, map_location='cpu')
        labels = labels_data['labels']
        side_effects = labels_data['side_effects']
        
        n_drugs, n_se = labels.shape
        log(f"\nOriginal: {n_drugs} drugs x {n_se} side effects", f)
        
        # Calculate positive rates
        pos_rates = labels.mean(dim=0)
        log(f"Original mean positive rate: {pos_rates.mean():.2%}", f)
        log(f"Original min: {pos_rates.min():.2%}, max: {pos_rates.max():.2%}", f)
        
        # === Strategy 1: Filter by positive rate (20-80%) ===
        lo, hi = 0.20, 0.80
        mask = (pos_rates >= lo) & (pos_rates <= hi)
        n_kept = mask.sum().item()
        log(f"\nSide effects with {lo:.0%}-{hi:.0%} positive rate: {n_kept}", f)
        
        # If too few, relax thresholds
        if n_kept < 30:
            lo, hi = 0.15, 0.85
            mask = (pos_rates >= lo) & (pos_rates <= hi)
            n_kept = mask.sum().item()
            log(f"Relaxed to {lo:.0%}-{hi:.0%}: {n_kept} side effects", f)
        
        if n_kept < 20:
            lo, hi = 0.10, 0.90
            mask = (pos_rates >= lo) & (pos_rates <= hi)
            n_kept = mask.sum().item()
            log(f"Relaxed to {lo:.0%}-{hi:.0%}: {n_kept} side effects", f)
        
        # Apply filter
        kept_indices = torch.where(mask)[0]
        filtered_labels = labels[:, kept_indices]
        filtered_se = [side_effects[i] for i in kept_indices.tolist()]
        
        new_pos_rates = filtered_labels.mean(dim=0)
        log(f"\nAfter filtering:", f)
        log(f"  Side effects kept: {len(filtered_se)}", f)
        log(f"  Mean positive rate: {new_pos_rates.mean():.2%}", f)
        log(f"  Min: {new_pos_rates.min():.2%}, Max: {new_pos_rates.max():.2%}", f)
        
        log(f"\nFiltered side effects:", f)
        for i, se in enumerate(filtered_se):
            rate = new_pos_rates[i].item()
            log(f"  {se}: {rate:.1%}", f)
        
        # Save filtered labels
        filtered_data = {
            'labels': filtered_labels,
            'drug_ids': labels_data.get('drug_ids', list(range(n_drugs))),
            'side_effects': filtered_se,
            'se_to_idx': {se: i for i, se in enumerate(filtered_se)},
            'original_indices': kept_indices.tolist(),
        }
        
        save_path = DATA_DIR / 'side_effect_labels_balanced.pt'
        torch.save(filtered_data, save_path)
        log(f"\nSaved filtered labels to {save_path}", f)
        
        return filtered_labels, filtered_se

def augment_graph(filtered_labels, filtered_se):
    """Augment graph with synthetic drug nodes using embedding noise."""
    with open(OUTPUT_FILE, 'a', encoding='utf-8') as f:
        log("\n" + "=" * 60, f)
        log("PHASE 2: AUGMENTING GRAPH WITH SYNTHETIC DRUGS", f)
        log("=" * 60, f)
        
        # Load graph
        graph_data = torch.load(DATA_DIR / 'graph_data.pt', weights_only=False, map_location='cpu')
        
        original_drug_x = graph_data['drug'].x.clone()
        n_original = original_drug_x.shape[0]
        emb_dim = original_drug_x.shape[1]
        
        log(f"\nOriginal graph: {n_original} drugs, embedding dim={emb_dim}", f)
        
        # === Augmentation: Add noisy copies ===
        n_augmented_per_drug = 2
        noise_std = 0.1
        
        aug_embeddings = []
        aug_labels = []
        
        for i in range(n_original):
            drug_emb = original_drug_x[i]
            drug_label = filtered_labels[i]
            
            for _ in range(n_augmented_per_drug):
                # Add Gaussian noise to embedding
                noise = torch.randn_like(drug_emb) * noise_std
                aug_emb = drug_emb + noise
                aug_embeddings.append(aug_emb)
                
                # For augmented drugs, flip some labels probabilistically
                # Side effects with high positive rate -> more likely to flip to 0
                # This creates more negatives for common side effects
                flip_prob = filtered_labels.mean(dim=0) * 0.3  # Higher base rate = more flips
                flip_mask = torch.bernoulli(flip_prob).bool()
                aug_label = drug_label.clone()
                aug_label[flip_mask] = 1.0 - aug_label[flip_mask]  # Flip
                aug_labels.append(aug_label)
        
        aug_embeddings = torch.stack(aug_embeddings)
        aug_labels = torch.stack(aug_labels)
        
        log(f"Generated {len(aug_embeddings)} augmented drugs", f)
        
        # Combine original + augmented
        combined_x = torch.cat([original_drug_x, aug_embeddings], dim=0)
        combined_labels = torch.cat([filtered_labels, aug_labels], dim=0)
        
        log(f"Combined: {combined_x.shape[0]} total drugs", f)
        
        # New positive rates after augmentation
        new_pos_rates = combined_labels.mean(dim=0)
        log(f"Post-augmentation mean positive rate: {new_pos_rates.mean():.2%}", f)
        log(f"Post-augmentation min: {new_pos_rates.min():.2%}, max: {new_pos_rates.max():.2%}", f)
        
        # Create augmented graph (copy original structure, update drug features)
        from torch_geometric.data import HeteroData
        aug_graph = HeteroData()
        
        aug_graph['drug'].x = combined_x
        aug_graph['drug'].num_nodes = combined_x.shape[0]
        
        # Copy protein data if it exists
        if 'protein' in graph_data.node_types:
            aug_graph['protein'].x = graph_data['protein'].x
            aug_graph['protein'].num_nodes = graph_data['protein'].num_nodes
        
        # Copy edge indices
        for edge_type in graph_data.edge_types:
            aug_graph[edge_type].edge_index = graph_data[edge_type].edge_index
        
        # Save
        graph_save_path = DATA_DIR / 'graph_data_augmented.pt'
        torch.save(aug_graph, graph_save_path)
        log(f"\nSaved augmented graph to {graph_save_path}", f)
        
        labels_save_path = DATA_DIR / 'side_effect_labels_augmented.pt'
        torch.save({
            'labels': combined_labels,
            'side_effects': filtered_se,
            'se_to_idx': {se: i for i, se in enumerate(filtered_se)},
            'n_original': n_original,
            'n_augmented': len(aug_embeddings),
        }, labels_save_path)
        log(f"Saved augmented labels to {labels_save_path}", f)
        
        log(f"\n{'='*60}", f)
        log("PHASE 1 & 2 COMPLETE", f)
        log(f"{'='*60}", f)

if __name__ == "__main__":
    filtered_labels, filtered_se = rebalance_labels()
    augment_graph(filtered_labels, filtered_se)
