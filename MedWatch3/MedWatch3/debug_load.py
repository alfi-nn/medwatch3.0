import sys
import os
from pathlib import Path

# Add root to path
BASE_DIR = Path(__file__).resolve().parent
sys.path.append(str(BASE_DIR))

from scripts.predict_new_drug import load_models, predict_side_effects

try:
    print("Testing model load...")
    load_models()
    print("Model load successful!")
    
    print("Testing prediction on Aspirin...")
    res = predict_side_effects("CC(=O)OC1=CC=CC=C1C(=O)O")
    print("Prediction successful!")
    print(f"Results keys: {list(res.keys())[:5]}...")
    
except Exception as e:
    print(f"FAILED: {e}")
    import traceback
    traceback.print_exc()
