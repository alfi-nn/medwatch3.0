
import torch
import pandas as pd
import sys
from pathlib import Path

# Paths
BASE_DIR = Path(r"c:\Users\kjalf\Downloads\project\MedWatch3\MedWatch3")
DATA_DIR = BASE_DIR / 'data' / 'processed'

def analyze_data():
    print(f"Analyzing data in {DATA_DIR}...")
    
    # Load labels
    try:
        labels_data = torch.load(DATA_DIR / 'side_effect_labels.pt', weights_only=False)
        labels = labels_data['labels']
        side_effects = labels_data['side_effects']
        
        print(f"Loaded {len(side_effects)} side effects for {labels.shape[0]} drugs.")
        
        # Calculate positive rates
        pos_counts = labels.sum(dim=0)
        pos_rates = pos_counts / labels.shape[0]
        
        print("\nTop 20 most common side effects:")
        sorted_indices = pos_rates.argsort(descending=True)
        for i in range(min(20, len(side_effects))):
            idx = sorted_indices[i]
            print(f"{side_effects[idx]}: {pos_rates[idx]:.2%} ({int(pos_counts[idx])} positives)")
            
        print("\nTop 20 LEAST common side effects:")
        for i in range(min(20, len(side_effects))):
            idx = sorted_indices[-(i+1)]
            print(f"{side_effects[idx]}: {pos_rates[idx]:.2%} ({int(pos_counts[idx])} positives)")
            
        # Distribution stats
        print("\nDistribution statistics:")
        print(f"Mean positive rate: {pos_rates.mean():.2%}")
        print(f"Median positive rate: {pos_rates.median():.2%}")
        print(f"Max positive rate: {pos_rates.max():.2%}")
        print(f"Min positive rate: {pos_rates.min():.2%}")
        
        # Check for zero columns
        zero_cols = (pos_counts == 0).sum().item()
        print(f"\nSide effects with ZERO positives: {zero_cols}")
        
    except Exception as e:
        print(f"Error loading labels: {e}")

if __name__ == "__main__":
    analyze_data()
