"""
Data Collection Pipeline for Colab - Saves directly to Google Drive
"""
import os
import sys
import gzip
import time
import logging
from pathlib import Path
from typing import Set, Dict, List, Optional
from collections import defaultdict

import pandas as pd
import pubchempy as pcp
import requests
from tqdm import tqdm

# Configuration - will be set from notebook
DRIVE_RAW = Path(os.environ.get('DRIVE_DATA_RAW', 'data/raw'))
DRIVE_PROCESSED = Path(os.environ.get('DRIVE_DATA_PROCESSED', 'data/processed'))

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# API settings
API_RATE_LIMIT = 0.5
TOP_N_DRUGS = 500

def download_file(url: str, output_path: Path) -> bool:
    """Download file with retry logic."""
    if output_path.exists():
        logger.info(f"File exists: {output_path.name}")
        return True
    
    headers = {'User-Agent': 'Mozilla/5.0'}
    try:
        logger.info(f"Downloading {url}...")
        response = requests.get(url, headers=headers, stream=True, timeout=300)
        response.raise_for_status()
        
        total = int(response.headers.get('content-length', 0))
        with open(output_path, 'wb') as f:
            with tqdm(total=total, unit='B', unit_scale=True) as pbar:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
                    pbar.update(len(chunk))
        logger.info(f"✅ Downloaded: {output_path.name}")
        return True
    except Exception as e:
        logger.error(f"Download failed: {e}")
        return False

def download_sider():
    """Download SIDER database files."""
    base_url = "http://sideeffects.embl.de/media/download/"
    files = [
        ("drug_names.tsv", "drug_names.tsv"),
        ("meddra_all_se.tsv.gz", "meddra_all_se.tsv.gz"),
    ]
    
    for remote, local in files:
        download_file(base_url + remote, DRIVE_RAW / local)

def download_stitch():
    """Download STITCH protein-chemical links."""
    url = "http://stitch.embl.de/download/protein_chemical.links.v5.0/9606.protein_chemical.links.v5.0.tsv.gz"
    download_file(url, DRIVE_RAW / "stitch_protein_chemical_links.tsv.gz")

def get_top_drugs(n: int = 500) -> pd.DataFrame:
    """Get top N drugs from SIDER by side effect frequency."""
    se_file = DRIVE_RAW / "meddra_all_se.tsv.gz"
    
    logger.info("Loading SIDER side effects...")
    df = pd.read_csv(se_file, sep='\t', compression='gzip', header=None,
                     names=['stitch_flat', 'stitch_stereo', 'umls_label', 'meddra_type', 
                            'umls_concept', 'side_effect'])
    
    # Count side effects per drug
    drug_counts = df['stitch_flat'].value_counts().head(n)
    top_drugs = drug_counts.index.tolist()
    
    # Convert STITCH ID to PubChem CID
    drug_data = []
    for stitch_id in top_drugs:
        # Remove CID prefix and leading zeros
        cid = stitch_id.replace('CID', '').replace('m', '').replace('s', '').lstrip('0')
        drug_data.append({
            'drug_id': f'CID{cid}',
            'stitch_id': stitch_id,
            'pubchem_cid': cid
        })
    
    return pd.DataFrame(drug_data)

def fetch_smiles(pubchem_cids: List[str]) -> Dict[str, str]:
    """Fetch SMILES from PubChem for given CIDs."""
    smiles_map = {}
    
    for cid in tqdm(pubchem_cids, desc="Fetching SMILES"):
        try:
            time.sleep(API_RATE_LIMIT)
            compounds = pcp.get_compounds(cid, 'cid')
            if compounds:
                smiles_map[cid] = compounds[0].isomeric_smiles or compounds[0].canonical_smiles
        except Exception as e:
            logger.warning(f"Failed to fetch SMILES for {cid}: {e}")
    
    return smiles_map

def collect_ddi_from_sider() -> pd.DataFrame:
    """Collect drug-drug interactions from SIDER co-occurrence."""
    se_file = DRIVE_RAW / "meddra_all_se.tsv.gz"
    
    logger.info("Building DDI edges from SIDER co-occurrence...")
    df = pd.read_csv(se_file, sep='\t', compression='gzip', header=None,
                     names=['stitch_flat', 'stitch_stereo', 'umls_label', 'meddra_type',
                            'umls_concept', 'side_effect'])
    
    # Group drugs by shared side effects
    se_to_drugs = defaultdict(set)
    for _, row in df.iterrows():
        se_to_drugs[row['side_effect']].add(row['stitch_flat'])
    
    # Create edges for drugs sharing side effects
    edges = set()
    for drugs in tqdm(se_to_drugs.values(), desc="Building DDI edges"):
        drugs = list(drugs)
        for i in range(len(drugs)):
            for j in range(i+1, min(i+10, len(drugs))):  # Limit pairs
                d1 = drugs[i].replace('CID', '').replace('m', '').replace('s', '').lstrip('0')
                d2 = drugs[j].replace('CID', '').replace('m', '').replace('s', '').lstrip('0')
                edges.add((f'CID{d1}', f'CID{d2}'))
    
    return pd.DataFrame(list(edges), columns=['source_drug_id', 'target_drug_id'])

def run_pipeline():
    """Run the complete data collection pipeline."""
    logger.info("=" * 60)
    logger.info("Starting Data Collection Pipeline")
    logger.info("=" * 60)
    
    # Step 1: Download raw data
    logger.info("\n📥 Step 1: Downloading raw data...")
    download_sider()
    download_stitch()
    
    # Step 2: Get top drugs
    logger.info("\n💊 Step 2: Identifying top drugs...")
    df_drugs = get_top_drugs(TOP_N_DRUGS)
    
    # Step 3: Fetch SMILES
    logger.info("\n🧪 Step 3: Fetching SMILES from PubChem...")
    smiles_map = fetch_smiles(df_drugs['pubchem_cid'].tolist())
    df_drugs['smiles'] = df_drugs['pubchem_cid'].map(smiles_map)
    df_drugs = df_drugs.dropna(subset=['smiles'])
    
    # Save drugs
    df_drugs.to_csv(DRIVE_PROCESSED / 'nodes_drugs.csv', index=False)
    logger.info(f"✅ Saved {len(df_drugs)} drugs to nodes_drugs.csv")
    
    # Step 4: Collect DDI edges
    logger.info("\n🔗 Step 4: Building DDI edges...")
    df_ddi = collect_ddi_from_sider()
    df_ddi.to_csv(DRIVE_PROCESSED / 'edges_ddi.csv', index=False)
    logger.info(f"✅ Saved {len(df_ddi)} DDI edges")
    
    # Step 5: Create empty protein/DTI files (STITCH parsing is slow)
    logger.info("\n🧬 Step 5: Creating placeholder protein files...")
    pd.DataFrame(columns=['protein_id', 'amino_acid_sequence']).to_csv(
        DRIVE_PROCESSED / 'nodes_proteins.csv', index=False)
    pd.DataFrame(columns=['drug_id', 'protein_id']).to_csv(
        DRIVE_PROCESSED / 'edges_dti.csv', index=False)
    
    logger.info("\n" + "=" * 60)
    logger.info("✅ Data Collection Complete!")
    logger.info("=" * 60)

if __name__ == "__main__":
    run_pipeline()
