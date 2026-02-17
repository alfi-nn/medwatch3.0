"""
Build heterogeneous graph with BERT embeddings.
"""
import os
import torch
import pandas as pd
import numpy as np
from pathlib import Path
from torch_geometric.data import HeteroData

DRIVE_PROCESSED = Path(os.environ.get('DRIVE_DATA_PROCESSED', 'data/processed'))

def build_graph():
    print("Building heterogeneous graph...")
    data = HeteroData()
    
    # Load data
    df_drugs = pd.read_csv(DRIVE_PROCESSED / 'nodes_drugs.csv')
    df_ddi = pd.read_csv(DRIVE_PROCESSED / 'edges_ddi.csv')
    
    drug_id_map = {id: i for i, id in enumerate(df_drugs['drug_id'])}
    
    # Load embeddings
    emb_path = DRIVE_PROCESSED / 'drug_embeddings.pt'
    if emb_path.exists():
        emb_data = torch.load(emb_path)
        embeddings = emb_data['embeddings']
        saved_ids = emb_data['drug_ids']
        
        # Align embeddings
        id_to_emb = {id: embeddings[i] for i, id in enumerate(saved_ids)}
        aligned = [id_to_emb.get(id, torch.zeros(embeddings.shape[1])) 
                   for id in df_drugs['drug_id']]
        data['drug'].x = torch.stack(aligned)
        print(f"✅ Using ChemBERTa embeddings: {data['drug'].x.shape}")
    else:
        print("⚠️ No embeddings found, using random features")
        data['drug'].x = torch.randn(len(df_drugs), 768)
    
    data['drug'].num_nodes = len(df_drugs)
    data['drug'].node_ids = df_drugs['drug_id'].tolist()
    
    # DDI edges
    src, dst = [], []
    for _, row in df_ddi.iterrows():
        s, d = row['source_drug_id'], row['target_drug_id']
        if s in drug_id_map and d in drug_id_map:
            src.append(drug_id_map[s])
            dst.append(drug_id_map[d])
    
    data['drug', 'interacts', 'drug'].edge_index = torch.tensor([src, dst], dtype=torch.long)
    print(f"✅ DDI edges: {len(src)}")
    
    # Save
    output_path = DRIVE_PROCESSED / 'graph_data.pt'
    torch.save(data, output_path)
    print(f"✅ Saved graph to {output_path}")
    print(data)

if __name__ == "__main__":
    build_graph()
