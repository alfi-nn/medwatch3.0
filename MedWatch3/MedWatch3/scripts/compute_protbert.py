"""Compute ProtBERT embeddings for proteins."""
import os
import torch
import pandas as pd
from pathlib import Path
from tqdm import tqdm

DRIVE_PROCESSED = Path(os.environ.get('DRIVE_DATA_PROCESSED', '/content/drive/MyDrive/MedWatch3/data/processed'))
device = 'cuda' if torch.cuda.is_available() else 'cpu'
print(f"Using device: {device}")

def main():
    df = pd.read_csv(DRIVE_PROCESSED / 'nodes_proteins.csv')
    print(f"Processing {len(df)} proteins...")
    
    from transformers import AutoTokenizer, AutoModel
    tokenizer = AutoTokenizer.from_pretrained("Rostlab/prot_bert_bfd")
    model = AutoModel.from_pretrained("Rostlab/prot_bert_bfd").to(device)
    model.eval()
    
    embeddings = []
    seq_col = 'amino_acid_sequence' if 'amino_acid_sequence' in df.columns else 'sequence'
    
    for seq in tqdm(df[seq_col], desc="Encoding proteins"):
        # ProtBERT expects spaces between amino acids
        seq_clean = str(seq)[:500] if pd.notna(seq) else 'M'
        seq_spaced = ' '.join(list(seq_clean))
        
        inputs = tokenizer(seq_spaced, return_tensors='pt', 
                          truncation=True, max_length=512, padding=True).to(device)
        with torch.no_grad():
            outputs = model(**inputs)
            emb = outputs.last_hidden_state[:, 0, :].cpu()
        embeddings.append(emb)
    
    embeddings = torch.cat(embeddings, dim=0)
    
    # Use protein_id column
    pid_col = 'protein_id' if 'protein_id' in df.columns else 'uniprot_id'
    
    torch.save({
        'embeddings': embeddings,
        'protein_ids': df[pid_col].tolist(),
        'model': 'Rostlab/prot_bert_bfd'
    }, DRIVE_PROCESSED / 'protein_embeddings.pt')
    
    print(f"✅ Saved protein embeddings: {embeddings.shape}")

if __name__ == "__main__":
    main()
