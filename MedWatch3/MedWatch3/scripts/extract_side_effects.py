"""Extract side effect labels from SIDER - FIXED matching."""
import gzip
import pandas as pd
import torch
from pathlib import Path
from collections import defaultdict

DRIVE_RAW = Path('/content/drive/MyDrive/MedWatch3/data/raw')
DRIVE_PROCESSED = Path('/content/drive/MyDrive/MedWatch3/data/processed')

def main():
    print("Loading SIDER side effects...")
    
    se_file = DRIVE_RAW / 'meddra_all_se.tsv.gz'
    df_se = pd.read_csv(se_file, sep='\t', compression='gzip', header=None,
                        names=['stitch_flat', 'stitch_stereo', 'umls_label', 
                               'meddra_type', 'umls_concept', 'side_effect'])
    
    print(f"Total side effect records: {len(df_se)}")
    
    # Load our drugs
    df_drugs = pd.read_csv(DRIVE_PROCESSED / 'nodes_drugs.csv')
    our_drug_ids = set(df_drugs['drug_id'].tolist())
    
    print(f"Our drugs: {len(our_drug_ids)}")
    print(f"Sample our IDs: {list(our_drug_ids)[:3]}")
    print(f"Sample SIDER IDs: {df_se['stitch_flat'].head(3).tolist()}")
    
    # Direct matching - IDs are already in same format!
    df_matched = df_se[df_se['stitch_flat'].isin(our_drug_ids)]
    print(f"Matched records: {len(df_matched)}")
    
    if len(df_matched) == 0:
        print("❌ No direct matches. Trying alternative...")
        # Maybe SIDER uses stitch_stereo column
        df_matched = df_se[df_se['stitch_stereo'].isin(our_drug_ids)]
        print(f"Matched via stitch_stereo: {len(df_matched)}")
    
    # Get top 100 side effects from matched data
    if len(df_matched) > 0:
        se_counts = df_matched['side_effect'].value_counts()
        top_side_effects = se_counts.head(100).index.tolist()
    else:
        # Fallback: use global top side effects
        top_side_effects = df_se['side_effect'].value_counts().head(100).index.tolist()
    
    print(f"Using top {len(top_side_effects)} side effects")
    se_to_idx = {se: i for i, se in enumerate(top_side_effects)}
    
    # Build drug -> side effects mapping
    drug_side_effects = defaultdict(set)
    for _, row in df_matched.iterrows():
        drug_id = row['stitch_flat']
        side_effect = row['side_effect']
        if side_effect in se_to_idx:
            drug_side_effects[drug_id].add(side_effect)
    
    print(f"Drugs with side effects: {len(drug_side_effects)}")
    
    # Create label matrix
    drug_ids = df_drugs['drug_id'].tolist()
    labels = torch.zeros(len(drug_ids), len(top_side_effects))
    
    for i, drug_id in enumerate(drug_ids):
        for se in drug_side_effects.get(drug_id, []):
            labels[i, se_to_idx[se]] = 1
    
    # Save
    torch.save({
        'labels': labels,
        'drug_ids': drug_ids,
        'side_effects': top_side_effects,
        'se_to_idx': se_to_idx,
    }, DRIVE_PROCESSED / 'side_effect_labels.pt')
    
    print(f"\n✅ Saved side effect labels: {labels.shape}")
    print(f"   Positive labels: {labels.sum().item():.0f}")
    print(f"\nTop 10 side effects: {top_side_effects[:10]}")

if __name__ == "__main__":
    main()
