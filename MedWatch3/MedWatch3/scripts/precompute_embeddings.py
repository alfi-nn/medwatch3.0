"""
Precompute ChemBERTa and ProtBERT embeddings.
Saves directly to Google Drive.
"""
import os
import torch
import pandas as pd
from pathlib import Path
from tqdm import tqdm

DRIVE_PROCESSED = Path(os.environ.get('DRIVE_DATA_PROCESSED', 'data/processed'))
device = 'cuda' if torch.cuda.is_available() else 'cpu'
print(f"Using device: {device}")

def compute_drug_embeddings():
    """Compute ChemBERTa embeddings for drugs."""
    output_path = DRIVE_PROCESSED / 'drug_embeddings.pt'
    
    if output_path.exists():
        print(f"✅ Drug embeddings already exist: {output_path}")
        return
    
    print("Loading ChemBERTa model...")
    from transformers import AutoTokenizer, AutoModel
    
    tokenizer = AutoTokenizer.from_pretrained("DeepChem/ChemBERTa-77M-MLM")
    model = AutoModel.from_pretrained("DeepChem/ChemBERTa-77M-MLM").to(device)
    model.eval()
    
    # Load drugs
    df = pd.read_csv(DRIVE_PROCESSED / 'nodes_drugs.csv')
    print(f"Processing {len(df)} drugs...")
    
    embeddings = []
    batch_size = 32
    
    for i in tqdm(range(0, len(df), batch_size)):
        batch = df['smiles'].iloc[i:i+batch_size].tolist()
        batch = [s if pd.notna(s) else 'C' for s in batch]
        
        inputs = tokenizer(batch, return_tensors='pt', padding=True, 
                          truncation=True, max_length=512).to(device)
        
        with torch.no_grad():
            outputs = model(**inputs)
            emb = outputs.last_hidden_state[:, 0, :].cpu()
        embeddings.append(emb)
    
    embeddings = torch.cat(embeddings, dim=0)
    
    torch.save({
        'embeddings': embeddings,
        'drug_ids': df['drug_id'].tolist(),
        'embedding_dim': embeddings.shape[1],
        'model': 'DeepChem/ChemBERTa-77M-MLM'
    }, output_path)
    
    print(f"✅ Saved drug embeddings: {embeddings.shape}")

if __name__ == "__main__":
    compute_drug_embeddings()
    print("\n✅ Embeddings saved to Drive!")
