# MedWatch3 — AI-Powered Adverse Drug Reaction Prediction

**Predicting Drug Side Effects Before They Happen.**

MedWatch3 is an advanced AI platform designed to predict adverse drug reactions (ADRs) and side effects of pharmaceutical compounds using only their molecular structure. By leveraging Heterogeneous Graph Transformers (HGT) and transformer-based molecular embeddings (ChemBERTa), MedWatch3 analyzes complex biological interactions to forecast potential risks early in the drug discovery process.

## 🚀 Key Features

*   **Deep Learning Prediction**: Uses an ensemble of HGT models and ChemBERTa embeddings to predict the probability of 76 different side effects.
*   **Graph Neural Network**: Models drugs, proteins, and their interactions (DDI, DTI) in a heterogeneous knowledge graph.
*   **Real-time Analysis**: Enter a SMILES string to get instant risk assessments.
*   **Interactive Dashboard**: A premium, glassmorphism-styled web interface for visualizing risk profiles and prediction confidence.
*   **Full-Stack Architecture**: Python FastAPI backend serving a modern HTML/JS frontend.

## 🏗️ System Architecture

The core of MedWatch3 is a **Heterogeneous Graph Transformer (HGT)** that learns from a biomedical knowledge graph containing:
*   **Drugs** (Nodes): Initialized with ChemBERTa-77M-MLM embeddings (768-dim).
*   **Proteins** (Nodes): Initialized with ProtBERT embeddings (1024-dim).
*   **Relationships** (Edges): Drug-Drug Interactions (SIDER) and Drug-Target Interactions (ChEMBL).

The model uses an **ensemble approach** (V3 Optimized + Balanced V2) with temperature scaling to provide calibrated, reliable probabilities.

## 📂 Project Structure

```
project/
├── MedWatch3/
│   └── MedWatch3/
│       ├── api.py                  # FastAPI Backend Server
│       ├── scripts/                # Data processing & training scripts
│       │   ├── predict_new_drug.py # CLI Prediction Script
│       │   ├── build_graph.py      # Graph construction
│       │   └── ...
│       ├── data/                   # Raw and processed datasets
│       ├── checkpoints/            # Trained model weights (.pt files)
│       └── requirements.txt        # Python dependencies
├── medfront/                       # Frontend Web Application
│   ├── index.html                  # Landing Page
│   ├── analysis.html               # Analysis Dashboard
│   ├── assets/                     # Images and Videos
│   └── ...
└── README.md                       # Project Documentation
```

## 🛠️ Installation & Setup

### Prerequisites
*   **Python 3.8+**
*   **pip** (Python package installer)

### 1. Clone the Repository
```bash
git clone <repository-url>
cd project
```

### 2. Backend Setup
Navigate to the backend directory and set up the Python environment.

```bash
cd MedWatch3/MedWatch3

# Create a virtual environment (Windows)
py -m venv venv
# Activate the environment
venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

*Note: For PyTorch and PyTorch Geometric, ensure you install the versions compatible with your CUDA version (if using GPU) or CPU.*

### 3. Frontend Setup
The frontend is built with vanilla HTML/CSS/JS and does not require a build step. You can simply open the files in a browser, but for full API integration, the backend must be running.

## 🚀 Usage

### Starting the Backend API
From the `MedWatch3/MedWatch3` directory (with venv activated):

```bash
uvicorn api:app --reload
```
The API will start at `http://127.0.0.1:8000`.
*   **Health Check**: `GET /`
*   **Predict Endpoint**: `POST /predict`
*   **Stats Endpoint**: `GET /stats`

### Launching the Web Interface
1.  Ensure the API server is running.
2.  Navigate to the `medfront` folder.
3.  Open `index.html` in your web browser.
4.  Enter a SMILES string (e.g., Aspirin: `CC(=O)OC1=CC=CC=C1C(=O)O`) and click **Analyze**.

### CLI Prediction (Optional)
You can also run predictions directly from the command line:

```bash
# In MedWatch3/MedWatch3/
python scripts/predict_new_drug.py --smiles "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
```

## 📊 Model Performance

| Metric | Score | Notes |
| :--- | :--- | :--- |
| **AUC-ROC** | **0.76** | Ensemble Model |
| **Side Effects** | 76 | Modeled Classes |
| **Drugs** | 500+ | Training Set (SIDER) |

The final model uses **Temperature Scaling (T=3.5)** to ensure predicted probabilities are realistic (e.g., 60% confidence strictly implies a 60% chance of the side effect), avoiding overconfident predictions common in deep learning models.

## 🧪 Tech Stack

*   **Frontend**: HTML5, CSS3 (Glassmorphism), JavaScript (Vanilla)
*   **Backend**: Python, FastAPI, Uvicorn
*   **AI/ML**: PyTorch, PyTorch Geometric, Hugging Face Transformers
*   **Data Sources**: SIDER, ChEMBL, PubChem, STITCH

## 📜 Credits

Developed as a comprehensive solution for AI-driven pharmacovigilance.
Based on research into Heterogeneous Graph Transformers for biomedical link prediction.
