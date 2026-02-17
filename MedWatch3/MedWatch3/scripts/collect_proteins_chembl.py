"""Fetch drug-protein interactions from ChEMBL."""
import os
import time
import requests
import pandas as pd
from pathlib import Path
from tqdm import tqdm

DRIVE_PROCESSED = Path(os.environ.get('DRIVE_DATA_PROCESSED', 'data/processed'))

def fetch_targets_for_drug(pubchem_cid):
    """Fetch protein targets from ChEMBL for a PubChem CID."""
    try:
        # First, get ChEMBL ID from PubChem CID
        url = f"https://www.ebi.ac.uk/unichem/rest/src_compound_id/{pubchem_cid}/22"
        resp = requests.get(url, timeout=10)
        if not resp.ok:
            return []
        
        data = resp.json()
        if not data:
            return []
        
        chembl_id = data[0].get('src_compound_id')
        if not chembl_id:
            return []
        
        # Get targets for this ChEMBL compound
        url = f"https://www.ebi.ac.uk/chembl/api/data/mechanism?molecule_chembl_id={chembl_id}&format=json"
        resp = requests.get(url, timeout=10)
        if not resp.ok:
            return []
        
        data = resp.json()
        targets = []
        for mech in data.get('mechanisms', []):
            target_id = mech.get('target_chembl_id')
            if target_id:
                targets.append(target_id)
        
        return targets
    except:
        return []

def fetch_protein_sequence(target_id):
    """Fetch protein sequence from ChEMBL target."""
    try:
        url = f"https://www.ebi.ac.uk/chembl/api/data/target/{target_id}?format=json"
        resp = requests.get(url, timeout=10)
        if not resp.ok:
            return None
        
        data = resp.json()
        components = data.get('target_components', [])
        if components:
            seq = components[0].get('sequence')
            return seq[:1000] if seq else None
    except:
        return None

def main():
    df_drugs = pd.read_csv(DRIVE_PROCESSED / 'nodes_drugs.csv')
    print(f"Processing {len(df_drugs)} drugs...")
    
    dti_edges = []
    proteins = {}
    
    # Process first 100 drugs to save time
    for _, row in tqdm(df_drugs.head(100).iterrows(), total=100, desc="Fetching targets"):
        time.sleep(0.3)
        cid = row['pubchem_cid']
        drug_id = row['drug_id']
        
        targets = fetch_targets_for_drug(cid)
        for target in targets[:3]:  # Max 3 targets per drug
            dti_edges.append({'drug_id': drug_id, 'protein_id': target})
            if target not in proteins:
                seq = fetch_protein_sequence(target)
                if seq:
                    proteins[target] = seq
    
    print(f"\n✅ Found {len(dti_edges)} DTI edges")
    print(f"✅ Got {len(proteins)} unique proteins")
    
    if proteins:
        df_proteins = pd.DataFrame([
            {'protein_id': k, 'sequence': v} for k, v in proteins.items()
        ])
        df_proteins.to_csv(DRIVE_PROCESSED / 'nodes_proteins.csv', index=False)
        
        df_dti = pd.DataFrame(dti_edges)
        df_dti.to_csv(DRIVE_PROCESSED / 'edges_dti.csv', index=False)
        print("✅ Saved to Drive!")
    else:
        print("❌ No protein data found")

if __name__ == "__main__":
    main()
