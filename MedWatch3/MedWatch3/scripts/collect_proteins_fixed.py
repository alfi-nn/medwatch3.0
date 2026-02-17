"""Fixed: Fetch drug-protein interactions from PubChem."""
import os
import time
import requests
import pandas as pd
from pathlib import Path
from tqdm import tqdm

DRIVE_PROCESSED = Path(os.environ.get('DRIVE_DATA_PROCESSED', 'data/processed'))

def get_gene_ids(cid):
    """Get gene IDs from PubChem - FIXED parsing."""
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/xrefs/GeneID/JSON"
        resp = requests.get(url, timeout=10)
        if resp.ok:
            data = resp.json()
            # Correct path to gene IDs
            info_list = data.get('InformationList', {}).get('Information', [])
            if info_list:
                genes = info_list[0].get('GeneID', [])
                return genes[:5]  # Max 5 genes per drug
    except Exception as e:
        pass
    return []

def main():
    df_drugs = pd.read_csv(DRIVE_PROCESSED / 'nodes_drugs.csv')
    print(f"Processing {len(df_drugs)} drugs...")
    
    dti_edges = []
    genes_found = set()
    drugs_with_targets = 0
    
    for _, row in tqdm(df_drugs.iterrows(), total=len(df_drugs), desc="Fetching gene targets"):
        time.sleep(0.25)  # Rate limit
        cid = str(row['pubchem_cid'])
        drug_id = row['drug_id']
        
        genes = get_gene_ids(cid)
        if genes:
            drugs_with_targets += 1
            for gene_id in genes:
                gene_str = f"GENE{gene_id}"
                dti_edges.append({'drug_id': drug_id, 'protein_id': gene_str})
                genes_found.add(gene_str)
    
    print(f"\n📊 Results:")
    print(f"   Drugs with targets: {drugs_with_targets}/{len(df_drugs)}")
    print(f"   Total DTI edges: {len(dti_edges)}")
    print(f"   Unique genes: {len(genes_found)}")
    
    if dti_edges:
        # Save DTI edges
        df_dti = pd.DataFrame(dti_edges)
        df_dti.to_csv(DRIVE_PROCESSED / 'edges_dti.csv', index=False)
        
        # Create protein nodes (using gene ID as placeholder sequence)
        df_proteins = pd.DataFrame([
            {'protein_id': g, 'sequence': 'M' + 'A' * 99}  # Placeholder 100 AA
            for g in genes_found
        ])
        df_proteins.to_csv(DRIVE_PROCESSED / 'nodes_proteins.csv', index=False)
        
        print(f"\n✅ Saved to Drive!")
        print(f"   - edges_dti.csv: {len(dti_edges)} edges")
        print(f"   - nodes_proteins.csv: {len(genes_found)} proteins")
    else:
        print("\n❌ No DTI data found")

if __name__ == "__main__":
    main()
