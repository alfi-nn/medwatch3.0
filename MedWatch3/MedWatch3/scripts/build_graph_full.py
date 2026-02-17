"""Build graph with drugs, proteins, DDI and DTI edges."""
import os
import torch
import pandas as pd
from pathlib import Path
from torch_geometric.data import HeteroData

DRIVE_PROCESSED = Path('/content/drive/MyDrive/MedWatch3/data/processed')

def build_graph():
    print("Building heterogeneous graph with DTI edges...")
    data = HeteroData()
    
    # Load data
    df_drugs = pd.read_csv(DRIVE_PROCESSED / 'nodes_drugs.csv')
    df_proteins = pd.read_csv(DRIVE_PROCESSED / 'nodes_proteins.csv')
    df_ddi = pd.read_csv(DRIVE_PROCESSED / 'edges_ddi.csv')
    df_dti = pd.read_csv(DRIVE_PROCESSED / 'edges_dti.csv')
    
    # Create ID mappings
    drug_id_map = {id: i for i, id in enumerate(df_drugs['drug_id'])}
    protein_id_map = {id: i for i, id in enumerate(df_proteins['protein_id'])}
    
    # Load drug embeddings
    drug_emb = torch.load(DRIVE_PROCESSED / 'drug_embeddings.pt', weights_only=False)
    drug_emb_dict = {id: drug_emb['embeddings'][i] for i, id in enumerate(drug_emb['drug_ids'])}
    drug_features = torch.stack([drug_emb_dict.get(id, torch.zeros(384)) for id in df_drugs['drug_id']])
    
    data['drug'].x = drug_features
    data['drug'].num_nodes = len(df_drugs)
    print(f"✅ Drugs: {data['drug'].x.shape}")
    
    # Load protein embeddings
    prot_emb = torch.load(DRIVE_PROCESSED / 'protein_embeddings.pt', weights_only=False)
    prot_emb_dict = {id: prot_emb['embeddings'][i] for i, id in enumerate(prot_emb['protein_ids'])}
    protein_features = torch.stack([prot_emb_dict.get(id, torch.zeros(1024)) for id in df_proteins['protein_id']])
    
    data['protein'].x = protein_features
    data['protein'].num_nodes = len(df_proteins)
    print(f"✅ Proteins: {data['protein'].x.shape}")
    
    # DDI edges
    src, dst = [], []
    for _, row in df_ddi.iterrows():
        s, d = row.iloc[0], row.iloc[1]
        if s in drug_id_map and d in drug_id_map:
            src.append(drug_id_map[s])
            dst.append(drug_id_map[d])
    data['drug', 'interacts', 'drug'].edge_index = torch.tensor([src, dst], dtype=torch.long)
    print(f"✅ DDI edges: {len(src)}")
    
    # DTI edges
    src, dst = [], []
    for _, row in df_dti.iterrows():
        drug_id, prot_id = row['drug_id'], row['protein_id']
        if drug_id in drug_id_map and prot_id in protein_id_map:
            src.append(drug_id_map[drug_id])
            dst.append(protein_id_map[prot_id])
    data['drug', 'targets', 'protein'].edge_index = torch.tensor([src, dst], dtype=torch.long)
    # Add reverse edges
    data['protein', 'targeted_by', 'drug'].edge_index = torch.tensor([dst, src], dtype=torch.long)
    print(f"✅ DTI edges: {len(src)}")
    
    # Save
    torch.save(data, DRIVE_PROCESSED / 'graph_data.pt')
    print(f"\n✅ Saved graph!")
    print(data)

if __name__ == "__main__":
    build_graph()
