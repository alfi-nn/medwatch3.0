"""
HGT ADR Prediction - Data Collection Pipeline
==============================================

A comprehensive, production-ready pipeline to fetch, clean, and align biological data
from public databases (SIDER, PubChem, UniProt, STITCH) to create heterogeneous graph
nodes and edges for Adverse Drug Reaction prediction using Heterogeneous Graph Transformers.

Author: B.Tech Final Year Project
Purpose: Data preprocessing for HGT-based ADR prediction model

Output Files:
1. nodes_drugs.csv - [drug_id, pubchem_cid, smiles, drug_name]
2. nodes_proteins.csv - [protein_id, uniprot_id, amino_acid_sequence]
3. edges_dti.csv - [drug_id, protein_id, binding_affinity_score]
4. edges_ddi.csv - [source_drug_id, target_drug_id, interaction_type]
"""

import os
import re
import time
import gzip
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
from urllib.request import urlretrieve
from collections import defaultdict

import pandas as pd
import pubchempy as pcp
import requests
from tqdm import tqdm
from bioservices import UniProt
from Bio import Entrez

# Set email for NCBI Entrez (required)
Entrez.email = "your_email@example.com"


def setup_logging():
    """Setup logging after directories are created."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler('data/logs/pipeline.log'),
            logging.StreamHandler()
        ]
    )


logger = logging.getLogger(__name__)


class Config:
    """Configuration class for managing paths, URLs, and parameters."""
    
    # Data directories
    BASE_DIR = Path("data")
    RAW_DIR = BASE_DIR / "raw"
    PROCESSED_DIR = BASE_DIR / "processed"
    LOGS_DIR = BASE_DIR / "logs"
    
    # SIDER dataset URLs (SIDER 4.1)
    SIDER_BASE_URL = "http://sideeffects.embl.de/media/download/"
    SIDER_DRUG_NAMES = "drug_names.tsv"
    SIDER_SIDE_EFFECTS = "meddra_all_se.tsv.gz"  # Gzipped file
    SIDER_DRUG_ATCS = "drug_atc.tsv"
    
    # STITCH database URLs (for drug-protein interactions)
    STITCH_BASE_URL = "http://stitch.embl.de/download/"
    STITCH_VERSION = "v5.0"
    STITCH_PROTEIN_LINKS = f"protein_chemical.links.v5.0/9606.protein_chemical.links.v5.0.tsv.gz"
    STITCH_ACTIONS = f"actions.v5.0/9606.actions.v5.0.tsv.gz"
    
    # API parameters
    TOP_N_DRUGS = 500  # Top N most common drugs from SIDER
    API_RETRY_MAX = 3
    API_RETRY_DELAY = 2  # seconds
    API_RATE_LIMIT_DELAY = 0.5  # seconds between requests
    
    # UniProt parameters
    UNIPROT_BATCH_SIZE = 100
    
    @classmethod
    def setup_directories(cls):
        """Create necessary directories if they don't exist."""
        for directory in [cls.RAW_DIR, cls.PROCESSED_DIR, cls.LOGS_DIR]:
            directory.mkdir(parents=True, exist_ok=True)


class SIDERDataLoader:
    """Handles downloading and parsing SIDER database files."""
    
    @staticmethod
    def download_file(url: str, output_path: Path) -> bool:
        """Download a file with error handling."""
        try:
            if output_path.exists():
                logger.info(f"File already exists: {output_path.name}")
                return True
            
            logger.info(f"Downloading {output_path.name}...")
            urlretrieve(url, output_path)
            logger.info(f"Downloaded: {output_path.name}")
            return True
        except Exception as e:
            logger.error(f"Failed to download {url}: {str(e)}")
            return False
    
    @staticmethod
    def download_sider_dataset() -> bool:
        """Download SIDER 4.1 dataset files."""
        files_to_download = {
            Config.SIDER_DRUG_NAMES: Config.SIDER_BASE_URL + "drug_names.tsv",
            Config.SIDER_SIDE_EFFECTS: Config.SIDER_BASE_URL + "meddra_all_se.tsv.gz",
            Config.SIDER_DRUG_ATCS: Config.SIDER_BASE_URL + "drug_atc.tsv"
        }
        
        success = True
        for filename, url in files_to_download.items():
            output_path = Config.RAW_DIR / filename
            if not SIDERDataLoader.download_file(url, output_path):
                success = False
        
        return success
    
    @staticmethod
    def load_drug_names() -> pd.DataFrame:
        """Load and parse SIDER drug names."""
        file_path = Config.RAW_DIR / Config.SIDER_DRUG_NAMES
        df = pd.read_csv(file_path, sep='\t', header=None, 
                        names=['stitch_flat_id', 'stitch_stereo_id', 'drug_name'])
        logger.info(f"Loaded {len(df)} drug names from SIDER")
        return df
    
    @staticmethod
    def load_side_effects() -> pd.DataFrame:
        """Load and parse SIDER side effects."""
        file_path = Config.RAW_DIR / Config.SIDER_SIDE_EFFECTS
        # Read gzipped file
        df = pd.read_csv(file_path, sep='\t', header=None, compression='gzip',
                        names=['stitch_flat_id', 'stitch_stereo_id', 'umls_cui', 
                               'meddra_type', 'umls_name', 'side_effect_name'])
        logger.info(f"Loaded {len(df)} side effect records")
        return df
    
    @staticmethod
    def get_top_n_drugs(n: int = 500) -> pd.DataFrame:
        """Get top N drugs by frequency in SIDER side effects database."""
        side_effects_df = SIDERDataLoader.load_side_effects()
        drug_names_df = SIDERDataLoader.load_drug_names()
        
        # Count drug occurrences
        drug_counts = side_effects_df['stitch_flat_id'].value_counts().head(n)
        top_drug_ids = drug_counts.index.tolist()
        
        # Filter drug names
        top_drugs = drug_names_df[drug_names_df['stitch_flat_id'].isin(top_drug_ids)].copy()
        top_drugs = top_drugs.drop_duplicates(subset='stitch_flat_id')
        
        logger.info(f"Selected top {len(top_drugs)} drugs from SIDER")
        return top_drugs


class STITCHDataLoader:
    """Handles downloading and parsing STITCH database for drug-protein interactions."""
    
    @staticmethod
    def download_stitch_files() -> bool:
        """Download STITCH protein-chemical links."""
        url = Config.STITCH_BASE_URL + Config.STITCH_PROTEIN_LINKS
        output_path = Config.RAW_DIR / "stitch_protein_chemical_links.tsv.gz"
        
        if output_path.exists():
            logger.info("STITCH file already exists")
            return True
        
        try:
            logger.info("Downloading STITCH protein-chemical links (this may take a while)...")
            urlretrieve(url, output_path)
            logger.info("STITCH download complete")
            return True
        except Exception as e:
            logger.error(f"Failed to download STITCH: {str(e)}")
            return False
    
    @staticmethod
    def load_protein_chemical_links(drug_ids: Set[str]) -> pd.DataFrame:
        """Load STITCH protein-chemical links for specific drugs using chunked reading."""
        file_path = Config.RAW_DIR / "stitch_protein_chemical_links.tsv.gz"
        
        if not file_path.exists():
            logger.warning("STITCH file not found, attempting download...")
            STITCHDataLoader.download_stitch_files()
        
        try:
            # Generate all possible STITCH ID formats for our drugs
            stitch_drug_ids = set()
            for drug_id in drug_ids:
                # Remove 'CID' prefix and get numeric part
                cid_numeric = drug_id.replace('CID', '')
                
                # Try multiple formats:
                # 1. CIDm + numeric (with leading zeros as-is)
                stitch_drug_ids.add(f"CIDm{cid_numeric}")
                # 2. CIDs + numeric (stereo version)
                stitch_drug_ids.add(f"CIDs{cid_numeric}")
                # 3. CIDm + numeric without leading zeros
                try:
                    stitch_drug_ids.add(f"CIDm{str(int(cid_numeric))}")
                    stitch_drug_ids.add(f"CIDs{str(int(cid_numeric))}")
                except:
                    pass
                # 4. CIDm + zero-padded to 9 digits
                stitch_drug_ids.add(f"CIDm{cid_numeric.zfill(9)}")
                stitch_drug_ids.add(f"CIDs{cid_numeric.zfill(9)}")
            
            logger.info(f"Generated {len(stitch_drug_ids)} possible STITCH ID variants")
            logger.info(f"Loading STITCH protein-chemical links (chunked reading)...")
            
            # Read file in chunks to avoid memory issues and filter efficiently
            matching_chunks = []
            chunk_size = 500000  # Process 500K rows at a time
            total_rows = 0
            
            for chunk in pd.read_csv(file_path, sep='\t', compression='gzip', 
                                    chunksize=chunk_size):
                total_rows += len(chunk)
                # Filter this chunk for our drug IDs
                matching = chunk[chunk['chemical'].isin(stitch_drug_ids)]
                
                if len(matching) > 0:
                    matching_chunks.append(matching)
                    logger.info(f"Found {len(matching)} matches in chunk (total processed: {total_rows:,} rows)")
            
            logger.info(f"Finished scanning {total_rows:,} total STITCH rows")
            
            if matching_chunks:
                df = pd.concat(matching_chunks, ignore_index=True)
                logger.info(f"Loaded {len(df)} drug-protein interactions from STITCH")
                return df
            else:
                logger.warning("No STITCH data matched our drug IDs")
                return pd.DataFrame(columns=['chemical', 'protein', 'combined_score'])
                
        except Exception as e:
            logger.error(f"Failed to load STITCH data: {str(e)}")
            return pd.DataFrame()


class IDMapper:
    """Handles mapping between different biological database identifiers."""
    
    @staticmethod
    def stitch_to_pubchem(stitch_id: str) -> Optional[str]:
        """Convert STITCH flat ID to PubChem CID."""
        # SIDER uses format CID000012345, need to convert to 12345
        match = re.search(r'CID[sm]?0*(\d+)', stitch_id)
        if match:
            return match.group(1)
        return None
    
    @staticmethod
    def stitch_protein_to_uniprot(stitch_protein_id: str) -> Optional[str]:
        """Convert STITCH protein ID to UniProt ID."""
        # STITCH format: 9606.ENSP00000000000
        # Extract Ensembl protein ID and map to UniProt
        parts = stitch_protein_id.split('.')
        if len(parts) == 2 and parts[0] == '9606':
            ensembl_id = parts[1]
            return ensembl_id
        return None


class DrugFeatureExtractor:
    """Extracts molecular features for drugs using PubChem API."""
    
    @staticmethod
    def fetch_drug_smiles_batch(pubchem_cids: List[str]) -> Dict[str, Dict]:
        """Fetch SMILES and drug info for multiple PubChem CIDs."""
        results = {}
        
        for cid in tqdm(pubchem_cids, desc="Fetching drug SMILES"):
            try:
                time.sleep(Config.API_RATE_LIMIT_DELAY)
                
                # Retry logic
                for attempt in range(Config.API_RETRY_MAX):
                    try:
                        compound = pcp.Compound.from_cid(int(cid))
                        results[cid] = {
                            'smiles': compound.canonical_smiles,
                            'iupac_name': getattr(compound, 'iupac_name', None),
                            'molecular_formula': compound.molecular_formula
                        }
                        break
                    except Exception as e:
                        if attempt < Config.API_RETRY_MAX - 1:
                            time.sleep(Config.API_RETRY_DELAY * (attempt + 1))
                        else:
                            logger.warning(f"Failed to fetch CID {cid} after {Config.API_RETRY_MAX} attempts: {str(e)}")
                            results[cid] = {'smiles': None, 'iupac_name': None, 'molecular_formula': None}
            
            except Exception as e:
                logger.error(f"Error processing CID {cid}: {str(e)}")
                results[cid] = {'smiles': None, 'iupac_name': None, 'molecular_formula': None}
        
        return results


class ProteinFeatureExtractor:
    """Extracts protein sequences using UniProt API."""
    
    def __init__(self):
        self.uniprot = UniProt()
    
    def fetch_protein_sequence(self, uniprot_id: str) -> Optional[str]:
        """Fetch amino acid sequence for a UniProt ID."""
        try:
            time.sleep(Config.API_RATE_LIMIT_DELAY)
            
            for attempt in range(Config.API_RETRY_MAX):
                try:
                    result = self.uniprot.search(uniprot_id, columns="sequence", limit=1)
                    if result:
                        lines = result.strip().split('\n')
                        if len(lines) > 1:
                            return lines[1]  # Second line contains sequence
                    return None
                except Exception as e:
                    if attempt < Config.API_RETRY_MAX - 1:
                        time.sleep(Config.API_RETRY_DELAY * (attempt + 1))
                    else:
                        raise e
        except Exception as e:
            logger.warning(f"Failed to fetch sequence for {uniprot_id}: {str(e)}")
            return None
    
    def ensembl_to_uniprot(self, ensembl_id: str) -> Optional[str]:
        """Map Ensembl protein ID to UniProt ID."""
        try:
            time.sleep(Config.API_RATE_LIMIT_DELAY)
            
            # Search UniProt for Ensembl cross-reference
            query = f"database:(type:ensembl {ensembl_id})"
            result = self.uniprot.search(query, columns="id", limit=1)
            
            if result:
                lines = result.strip().split('\n')
                if len(lines) > 1:
                    return lines[1].split()[0]
            return None
        except Exception as e:
            logger.warning(f"Failed to map {ensembl_id} to UniProt: {str(e)}")
            return None


class InteractionCollector:
    """Collects drug-drug and drug-target interactions."""
    
    @staticmethod
    def collect_dti_from_stitch(drug_pubchem_map: Dict[str, str]) -> pd.DataFrame:
        """Collect Drug-Target Interactions from STITCH."""
        # Download STITCH data
        STITCHDataLoader.download_stitch_files()
        
        # Load STITCH data - pass drug IDs (with CID prefix), not just numeric CIDs
        drug_ids = set(drug_pubchem_map.keys())  # These are like "CID100000159"
        stitch_df = STITCHDataLoader.load_protein_chemical_links(drug_ids)
        
        if stitch_df.empty:
            logger.warning("No STITCH data available")
            return pd.DataFrame(columns=['drug_id', 'protein_id', 'binding_affinity_score'])
        
        # Process interactions
        dti_edges = []
        
        for _, row in tqdm(stitch_df.iterrows(), total=len(stitch_df), desc="Processing DTI edges"):
            chemical_id = row['chemical']
            protein_id = row['protein']
            score = row.get('combined_score', 500)  # STITCH confidence score
            
            # Convert STITCH IDs
            pubchem_cid = IDMapper.stitch_to_pubchem(chemical_id)
            ensembl_id = IDMapper.stitch_protein_to_uniprot(protein_id)
            
            if pubchem_cid and ensembl_id:
                # Find original drug ID
                drug_id = None
                for did, cid in drug_pubchem_map.items():
                    if cid == pubchem_cid:
                        drug_id = did
                        break
                
                if drug_id:
                    dti_edges.append({
                        'drug_id': drug_id,
                        'protein_id': ensembl_id,
                        'binding_affinity_score': score
                    })
        
        dti_df = pd.DataFrame(dti_edges)
        logger.info(f"Collected {len(dti_df)} DTI edges")
        return dti_df
    
    @staticmethod
    def collect_ddi_from_sider() -> pd.DataFrame:
        """Collect Drug-Drug Interactions from SIDER co-occurrence."""
        # Load side effects
        side_effects_df = SIDERDataLoader.load_side_effects()
        
        # Find drugs that share side effects (simple DDI proxy)
        ddi_edges = []
        
        # Group by side effect
        grouped = side_effects_df.groupby('umls_cui')['stitch_flat_id'].apply(list)
        
        logger.info("Generating DDI edges from shared side effects...")
        for side_effect, drugs in tqdm(grouped.items(), desc="Processing DDI edges"):
            if len(drugs) > 1:
                # Create edges between drugs sharing side effects
                drugs = list(set(drugs))
                for i in range(len(drugs)):
                    for j in range(i+1, len(drugs)):
                        ddi_edges.append({
                            'source_drug_id': drugs[i],
                            'target_drug_id': drugs[j],
                            'interaction_type': 'shared_side_effect'
                        })
        
        ddi_df = pd.DataFrame(ddi_edges)
        ddi_df = ddi_df.drop_duplicates(subset=['source_drug_id', 'target_drug_id'])
        
        logger.info(f"Collected {len(ddi_df)} DDI edges")
        return ddi_df


class GraphDataExporter:
    """Exports processed data to CSV files for PyTorch Geometric."""
    
    @staticmethod
    def export_nodes_drugs(drugs_df: pd.DataFrame, smiles_dict: Dict) -> pd.DataFrame:
        """Export nodes_drugs.csv."""
        output_path = Config.PROCESSED_DIR / "nodes_drugs.csv"
        
        # Create dataframe
        nodes = []
        for _, row in drugs_df.iterrows():
            drug_id = row['stitch_flat_id']
            pubchem_cid = IDMapper.stitch_to_pubchem(drug_id)
            
            if pubchem_cid and pubchem_cid in smiles_dict:
                nodes.append({
                    'drug_id': drug_id,
                    'pubchem_cid': pubchem_cid,
                    'smiles': smiles_dict[pubchem_cid].get('smiles'),
                    'drug_name': row['drug_name']
                })
        
        nodes_df = pd.DataFrame(nodes)
        nodes_df.to_csv(output_path, index=False)
        
        logger.info(f"Exported {len(nodes_df)} drug nodes to {output_path}")
        return nodes_df
    
    @staticmethod
    def export_nodes_proteins(protein_sequences: Dict[str, str]) -> pd.DataFrame:
        """Export nodes_proteins.csv."""
        output_path = Config.PROCESSED_DIR / "nodes_proteins.csv"
        
        nodes = []
        for protein_id, sequence in protein_sequences.items():
            if sequence:
                nodes.append({
                    'protein_id': protein_id,
                    'uniprot_id': protein_id,  # Will be updated with actual UniProt IDs
                    'amino_acid_sequence': sequence
                })
        
        nodes_df = pd.DataFrame(nodes)
        nodes_df.to_csv(output_path, index=False)
        
        logger.info(f"Exported {len(nodes_df)} protein nodes to {output_path}")
        return nodes_df
    
    @staticmethod
    def export_edges_dti(dti_df: pd.DataFrame) -> pd.DataFrame:
        """Export edges_dti.csv."""
        output_path = Config.PROCESSED_DIR / "edges_dti.csv"
        dti_df.to_csv(output_path, index=False)
        
        logger.info(f"Exported {len(dti_df)} DTI edges to {output_path}")
        return dti_df
    
    @staticmethod
    def export_edges_ddi(ddi_df: pd.DataFrame) -> pd.DataFrame:
        """Export edges_ddi.csv."""
        output_path = Config.PROCESSED_DIR / "edges_ddi.csv"
        ddi_df.to_csv(output_path, index=False)
        
        logger.info(f"Exported {len(ddi_df)} DDI edges to {output_path}")
        return ddi_df


def main():
    """Main pipeline orchestration."""
    print("=" * 70)
    print("HGT ADR Prediction - Data Collection Pipeline")
    print("=" * 70)
    print()
    
    # Setup
    Config.setup_directories()
    setup_logging()
    logger.info("Directories created and logging initialized")
    
    # Step 1: Download SIDER
    logger.info("\nStep 1/8: Downloading SIDER dataset...")
    if not SIDERDataLoader.download_sider_dataset():
        logger.error("Failed to download SIDER dataset. Exiting.")
        return
    
    # Step 2: Get top N drugs
    logger.info(f"\nStep 2/8: Extracting top {Config.TOP_N_DRUGS} drugs from SIDER...")
    top_drugs_df = SIDERDataLoader.get_top_n_drugs(Config.TOP_N_DRUGS)
    
    # Step 3: Map to PubChem and fetch SMILES
    logger.info("\nStep 3/8: Mapping to PubChem CIDs...")
    drug_pubchem_map = {}
    pubchem_cids = []
    
    for _, row in top_drugs_df.iterrows():
        stitch_id = row['stitch_flat_id']
        pubchem_cid = IDMapper.stitch_to_pubchem(stitch_id)
        if pubchem_cid:
            drug_pubchem_map[stitch_id] = pubchem_cid
            pubchem_cids.append(pubchem_cid)
    
    logger.info(f"Mapped {len(pubchem_cids)} drugs to PubChem CIDs")
    
    # Step 4: Fetch drug features
    logger.info("\nStep 4/8: Fetching drug SMILES from PubChem...")
    smiles_dict = DrugFeatureExtractor.fetch_drug_smiles_batch(pubchem_cids)
    
    # Step 5: Collect DTI edges
    logger.info("\nStep 5/8: Collecting Drug-Target Interactions...")
    dti_df = InteractionCollector.collect_dti_from_stitch(drug_pubchem_map)
    
    # Step 6: Fetch protein sequences
    logger.info("\nStep 6/8: Fetching protein sequences from UniProt...")
    protein_extractor = ProteinFeatureExtractor()
    protein_sequences = {}
    
    unique_proteins = dti_df['protein_id'].unique() if not dti_df.empty else []
    
    for ensembl_id in tqdm(unique_proteins[:1000], desc="Mapping proteins to UniProt"):  # Limit to avoid timeout
        uniprot_id = protein_extractor.ensembl_to_uniprot(ensembl_id)
        if uniprot_id:
            sequence = protein_extractor.fetch_protein_sequence(uniprot_id)
            if sequence:
                protein_sequences[uniprot_id] = sequence
    
    # Step 7: Collect DDI edges
    logger.info("\nStep 7/8: Collecting Drug-Drug Interactions...")
    ddi_df = InteractionCollector.collect_ddi_from_sider()
    
    # Filter DDI to only include drugs in our set
    drug_ids_set = set(drug_pubchem_map.keys())
    ddi_df = ddi_df[
        ddi_df['source_drug_id'].isin(drug_ids_set) & 
        ddi_df['target_drug_id'].isin(drug_ids_set)
    ]
    
    # Step 8: Export to CSV
    logger.info("\nStep 8/8: Exporting data to CSV files...")
    drugs_nodes_df = GraphDataExporter.export_nodes_drugs(top_drugs_df, smiles_dict)
    proteins_nodes_df = GraphDataExporter.export_nodes_proteins(protein_sequences)
    dti_edges_df = GraphDataExporter.export_edges_dti(dti_df)
    ddi_edges_df = GraphDataExporter.export_edges_ddi(ddi_df)
    
    # Summary
    print("\n" + "=" * 70)
    print("PIPELINE COMPLETE - Summary Statistics")
    print("=" * 70)
    print(f"Drug Nodes:       {len(drugs_nodes_df):,}")
    print(f"Protein Nodes:    {len(proteins_nodes_df):,}")
    print(f"DTI Edges:        {len(dti_edges_df):,}")
    print(f"DDI Edges:        {len(ddi_edges_df):,}")
    print()
    print(f"Output files saved to: {Config.PROCESSED_DIR.absolute()}")
    print("=" * 70)
    
    logger.info("Pipeline execution completed successfully!")


if __name__ == "__main__":
    main()
