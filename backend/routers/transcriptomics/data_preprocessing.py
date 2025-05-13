from fastapi import APIRouter, HTTPException
from fastapi.responses import JSONResponse
import os
import scanpy as sc
# import adata_preprocessor as ap
import hdf5plugin
from routers.transcriptomics.adata_preprocessor import adata_preprocessor as ap
router = APIRouter()

# Paths to directories

RAW_DATA_DIR = "routers/transcriptomics/data/1.H5AD_standalone/GSE161872"
PROCESSED_DATA_DIR = "routers/transcriptomics/data/2.H5AD_filtered/GSE161872"

@router.get("/preprocessing/summary")
def get_summary():
    try:
        raw_files = [f for f in os.listdir(RAW_DATA_DIR) if f.endswith(".h5ad")]
        processed_files = [f for f in os.listdir(PROCESSED_DATA_DIR) if f.endswith(".h5ad")]

        if not raw_files:
            raise HTTPException(status_code=404, detail="No raw .h5ad files found in the specified directory.")

        if not processed_files:
            for file in raw_files:
                raw_path = os.path.join(RAW_DATA_DIR, file)
                processed_path = os.path.join(PROCESSED_DATA_DIR, f"filtered_{file}")
                adata_filtered = ap.adata_preprocessor(raw_path, n_genes_min=400, n_genes_max=10000)
                adata_filtered.write_h5ad(processed_path, compression=hdf5plugin.FILTERS["zstd"])

            processed_files = [f for f in os.listdir(PROCESSED_DATA_DIR) if f.endswith(".h5ad")]

        summaries = []
        for raw_file, processed_file in zip(raw_files, processed_files):
            raw_path = os.path.join(RAW_DATA_DIR, raw_file)
            processed_path = os.path.join(PROCESSED_DATA_DIR, processed_file)

            # Load raw and processed data
            raw_data = sc.read_h5ad(raw_path)
            processed_data = sc.read_h5ad(processed_path)

            # Append summary for each pair
            summaries.append({
                "file_pair": {
                    "raw": raw_file,
                    "processed": processed_file
                },
                "raw_summary": {
                    "n_cells": raw_data.n_obs,
                    "n_genes": raw_data.n_vars,
                    "sparsity": 1.0 - (raw_data.X.nnz / (raw_data.shape[0] * raw_data.shape[1])) if hasattr(raw_data.X, 'nnz') else None
                },
                "processed_summary": {
                    "n_cells": processed_data.n_obs,
                    "n_genes": processed_data.n_vars,
                    "sparsity": 1.0 - (processed_data.X.nnz / (processed_data.shape[0] * processed_data.shape[1])) if hasattr(processed_data.X, 'nnz') else None
                }
            })

        return JSONResponse(content={"summaries": summaries})

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error generating summary: {str(e)}")

@router.get("/preprocessing/sample")
def get_sample(rows: int = 100, cols: int = 100):
    try:
        raw_files = [f for f in os.listdir(RAW_DATA_DIR) if f.endswith(".h5ad")]
        processed_files = [f for f in os.listdir(PROCESSED_DATA_DIR) if f.endswith(".h5ad")]

        if not raw_files:
            raise HTTPException(status_code=404, detail="No raw .h5ad files found in the specified directory.")

        if not processed_files:
            for file in raw_files:
                raw_path = os.path.join(RAW_DATA_DIR, file)
                processed_path = os.path.join(PROCESSED_DATA_DIR, f"filtered_{file}")
                adata_filtered = ap.adata_preprocessor(raw_path, n_genes_min=400, n_genes_max=10000)
                adata_filtered.write_h5ad(processed_path, compression=hdf5plugin.FILTERS["zstd"])

            processed_files = [f for f in os.listdir(PROCESSED_DATA_DIR) if f.endswith(".h5ad")]

        samples = []
        for raw_file, processed_file in zip(raw_files, processed_files):
            raw_path = os.path.join(RAW_DATA_DIR, raw_file)
            processed_path = os.path.join(PROCESSED_DATA_DIR, processed_file)

            # Load raw and processed data
            raw_data = sc.read_h5ad(raw_path)
            processed_data = sc.read_h5ad(processed_path)

            # Get a small sample of the data
            raw_sample = raw_data[:rows, :cols].to_df() if raw_data.shape[0] >= rows and raw_data.shape[1] >= cols else raw_data.to_df()
            processed_sample = processed_data[:rows, :cols].to_df() if processed_data.shape[0] >= rows and processed_data.shape[1] >= cols else processed_data.to_df()

            samples.append({
                "file_pair": {
                    "raw": raw_file,
                    "processed": processed_file
                },
                "raw_sample": raw_sample.to_dict(orient="split"),
                "processed_sample": processed_sample.to_dict(orient="split")
            })

        return JSONResponse(content={"samples": samples})

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error generating sample: {str(e)}")
