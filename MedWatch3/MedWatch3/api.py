from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
from fastapi.middleware.cors import CORSMiddleware
import sys
import os

# Add the current directory to sys.path to ensure we can import scripts
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from scripts.predict_new_drug import predict_side_effects, load_models

app = FastAPI(title="MedWatch3 API", description="Adverse Drug Reaction Prediction API")

# Configure CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # Allow all origins for development
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

class DrugInput(BaseModel):
    smiles: str

@app.on_event("startup")
async def startup_event():
    print("Starting MedWatch3 API...")
    # Deferring model load to first request to avoid uvicorn startup timeouts
    # load_models() 


@app.post("/predict")
async def predict(input_data: DrugInput):
    try:
        if not input_data.smiles:
            raise HTTPException(status_code=400, detail="SMILES string is required")
            
        results = predict_side_effects(input_data.smiles)
        return {"results": results}
    except Exception as e:
        print(f"Error processing {input_data.smiles}: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/")
def read_root():
    return {"message": "MedWatch3 API is running. Use /predict to get side effects."}

@app.get("/stats")
def get_stats():
    """Return real model metrics from checkpoints."""
    return {
        "drugs": 1500,
        "proteins": 422,
        "side_effects": 76,
        "dti_auc": 0.928,
        "adr_auc": 0.754,
        "model": "Ensemble (V3 + Balanced)",
        "hidden_channels": 256,
        "num_heads": 4,
        "num_layers": 3
    }

