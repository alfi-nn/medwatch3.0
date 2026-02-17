"""Parse STITCH and fetch protein sequences from UniProt."""
import os
import gzip
import time
import pandas as pd
from pathlib import Path
from tqdm import tqdm

DRIVE_RAW = Path(os.environ.get('DRIVE_DATA_RAW', 'data/raw'))
DRIVE_PROCESSED = Path(os.environ.get('DRIVE_DATA_PROCESSED', 'data/processed'))

def parse_stitch():
    """Parse STITCH protein-chemical links for our drugs."""
    print("📦 Loading existing drugs...")
    df_drugs = pd.read_csv(DRIVE_PROCESSED / 'nodes_drugs.csv')
    drug_ids = set(df_drugs['drug_id'].tolist())
    
    # Build STITCH ID variants
    stitch_ids = set()
    for drug_id in drug_ids:
        cid = drug_id.replace('CID', '')
        stitch_ids.add(f"CIDm{cid.zfill(9)}")
        stitch_ids.add(f"CIDs{cid.zfill(9)}")
    
    print(f"🔍 Scanning STITCH for {len(stitch_ids)} drug variants...")
    stitch_file = DRIVE_RAW / 'stitch_protein_chemical_links.tsv.gz'
    
    matches = []
    with gzip.open(stitch_file, 'rt') as f:
        header = f.readline()
        for line in tqdm(f, desc="Scanning STITCH"):
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                chemical, protein, score = parts[0], parts[1], parts[2]
                if chemical in stitch_ids and int(score) >= 700:
                    matches.append({
                        'stitch_id': chemical,
                        'protein_id': protein.replace('9606.', ''),
                        'score': int(score)
                    })
    
    print(f"✅ Found {len(matches)} drug-protein interactions")
    return pd.DataFrame(matches)

def fetch_sequences(protein_ids, max_proteins=200):
    """Fetch protein sequences from UniProt."""
    import requests
    
    unique_proteins = list(set(protein_ids))[:max_proteins]
    print(f"🧬 Fetching sequences for {len(unique_proteins)} proteins...")
    
    proteins = []
    for pid in tqdm(unique_proteins, desc="Fetching sequences"):
        try:
            time.sleep(0.5)
            # Convert Ensembl to UniProt
            url = f"https://rest.uniprot.org/uniprotkb/search?query={pid}&format=fasta&size=1"
            resp = requests.get(url, timeout=10)
            if resp.ok and resp.text.strip():
                lines = resp.text.strip().split('\n')
                seq = ''.join(lines[1:]) if len(lines) > 1 else None
                if seq and len(seq) > 50:
                    proteins.append({'protein_id': pid, 'sequence': seq[:1000]})
        except Exception as e:
            pass
    
    print(f"✅ Got sequences for {len(proteins)} proteins")
    return pd.DataFrame(proteins)

def main():
    # Step 1: Parse STITCH
    df_dti = parse_stitch()
    
    if len(df_dti) == 0:
        print("❌ No DTI interactions found!")
        return
    
    # Map STITCH IDs back to our drug IDs
    df_drugs = pd.read_csv(DRIVE_PROCESSED / 'nodes_drugs.csv')
    stitch_to_drug = {}
    for _, row in df_drugs.iterrows():
        cid = row['drug_id'].replace('CID', '')
        stitch_to_drug[f"CIDm{cid.zfill(9)}"] = row['drug_id']
        stitch_to_drug[f"CIDs{cid.zfill(9)}"] = row['drug_id']
    
    df_dti['drug_id'] = df_dti['stitch_id'].map(stitch_to_drug)
    df_dti = df_dti.dropna(subset=['drug_id'])
    
    # Step 2: Fetch sequences
    df_proteins = fetch_sequences(df_dti['protein_id'].tolist())
    
    # Step 3: Save
    df_proteins.to_csv(DRIVE_PROCESSED / 'nodes_proteins.csv', index=False)
    print(f"✅ Saved {len(df_proteins)} proteins")
    
    # Save DTI edges (only for proteins we have)
    valid_proteins = set(df_proteins['protein_id'])
    df_dti = df_dti[df_dti['protein_id'].isin(valid_proteins)]
    df_dti[['drug_id', 'protein_id']].to_csv(DRIVE_PROCESSED / 'edges_dti.csv', index=False)
    print(f"✅ Saved {len(df_dti)} DTI edges")

if __name__ == "__main__":
    main()
