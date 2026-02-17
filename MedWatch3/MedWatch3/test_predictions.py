
import sys
import torch
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent))

from scripts.predict_new_drug import predict_side_effects, print_predictions

def main():
    examples = {
        "Aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "Ibuprofen": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
    }
    
    with open('predictions_balanced_v2.txt', 'w', encoding='utf-8') as f:
        sys.stdout = f
        
        for name, smiles in examples.items():
            print(f"\n{'='*60}")
            print(f"Predicting side effects for {name} ({smiles})...")
            try:
                results = predict_side_effects(smiles)
                print_predictions(smiles, results)
            except Exception as e:
                print(f"Error predicting for {name}: {e}")
                import traceback
                traceback.print_exc()

if __name__ == "__main__":
    main()
