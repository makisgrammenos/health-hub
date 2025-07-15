from fastapi import APIRouter, HTTPException
from fastapi.responses import JSONResponse
import pandas as pd
from pycaret.clustering import setup, create_model, assign_model, plot_model, pull
import os

router = APIRouter()

@router.get("/clustering")
def run_clustering_demo():
    try:
        # Load dataset
        file_path = "routers/clinical/uci_hf_df.csv"
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"Dataset not found at {file_path}")
        
        df = pd.read_csv(file_path)
        
        # Set up PyCaret Clustering Experiment
        s = setup(df, session_id=123, ignore_features=["tag"], verbose=False)
        
        # Create and assign clusters using KMeans
        kmeans_model = create_model('kmeans', num_clusters=2)
        cluster_assignments = assign_model(kmeans_model)
        
        # Extract metrics
        metrics = pull()
        
        # Return cluster assignments and metrics
        return JSONResponse(content={
            "cluster_assignments": cluster_assignments.to_dict(orient="records"),
            "metrics": metrics.to_dict(orient="records")
        })

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error running clustering demo: {str(e)}")
