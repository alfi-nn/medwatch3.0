"""Fetch drug-protein interactions directly from PubChem."""
import os
import time
import requests
import pandas as pd
from pathlib import Path
from tqdm import tqdm

DRIVE_PROCESSED = Path(os.environ.get('DRIVE_DATA_PROCESSED', 'data/processed'))

def get_targets_from_pubchem(cid):
    """Get protein targets from PubChem for a compound."""
    try:
        # PubChem compound-gene link
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/xrefs/GeneID/JSON"
        resp = requests.get(url, timeout=10)
        if resp.ok:
            data = resp.json()
            genes = data.get('InformationList', {}).get('Information', [{}])[0].get('GeneID', [])
            return genes[:5]  # Max 5 genes per drug
    except:
        pass
    return []

def get_protein_sequence(gene_id):
    """Get protein sequence from NCBI for a gene ID."""
    try:
        # Get protein accession from gene
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id={gene_id}&rettype=fasta&retmode=text"
        resp = requests.get(url, timeout=10)
        if resp.ok and 'ORIGIN' in resp.text:
            # Extract sequence (simplified)
            lines = resp.text.split('\n')
            seq = ''.join([l for l in lines if not l.startswith('>')])[:500]
            return seq if len(seq) > 50 else None
    except:
        pass
    return None

def main():
    df_drugs = pd.read_csv(DRIVE_PROCESSED / 'nodes_drugs.csv')
    print(f"Processing {len(df_drugs)} drugs from PubChem...")
    
    dti_edges = []
    proteins = {}
    
    for _, row in tqdm(df_drugs.head(150).iterrows(), total=150, desc="Fetching from PubChem"):
        time.sleep(0.25)
        cid = row['pubchem_cid']
        drug_id = row['drug_id']
        
        genes = get_targets_from_pubchem(cid)
        for gene_id in genes:
            gene_str = f"GENE{gene_id}"
            dti_edges.append({'drug_id': drug_id, 'protein_id': gene_str})
            proteins[gene_str] = str(gene_id)  # Store gene ID for now
    
    print(f"\n✅ Found {len(dti_edges)} DTI edges")
    print(f"✅ Got {len(proteins)} unique genes/proteins")
    
    if dti_edges:
        # Save edges
        df_dti = pd.DataFrame(dti_edges)
        df_dti.to_csv(DRIVE_PROCESSED / 'edges_dti.csv', index=False)
        
        # For proteins, we'll use gene IDs (can fetch sequences later)
        df_proteins = pd.DataFrame([
            {'protein_id': k, 'sequence': 'PLACEHOLDER' * 50} for k in proteins.keys()
        ])
        df_proteins.to_csv(DRIVE_PROCESSED / 'nodes_proteins.csv', index=False)
        
        print("✅ Saved DTI edges and protein placeholders!")
        print("   (Using gene IDs as protein IDs)")
    else:
        print("❌ No data found")

if __name__ == "__main__":
    main()
