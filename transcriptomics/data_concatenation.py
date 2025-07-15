from fastapi import FastAPI, APIRouter
import scanpy as sc
import hdf5plugin
import anndata
import pandas as pd
import os

router = APIRouter()

VERSION = "1"
FILTERED_DIR = "routers/transcriptomics/data/2.H5AD_filtered/GSE161872"
OUTPUT_DIR = "routers/transcriptomics/data/3.H5AD_concatenated"
MERGED_NAME = "GSE161872"
H5AD_CONCAT = f"{OUTPUT_DIR}/{MERGED_NAME}_v{VERSION}.h5ad"

def get_h5ads(directory):
    """Fetch all H5AD files from the directory."""
    return [file for file in os.listdir(directory) if file.endswith(".h5ad")]

@router.post("/merge_h5ads")
def merge_h5ads():
    try:
        # Step 1: Retrieve the list of H5AD files
        adata_name_lst = get_h5ads(FILTERED_DIR)
        h5ad_files = [os.path.join(FILTERED_DIR, file) for file in adata_name_lst]

        # Step 2: Read and merge the H5AD files
        adatas = [sc.read_h5ad(file) for file in h5ad_files]
        adata_merged = anndata.AnnData.concatenate(*adatas, batch_key="batch_num", join="inner")

        # Step 3: Keep only 'gene_ids' in `adata_merged.var`
        adata_merged.var = adata_merged.var[["gene_ids"]]

        # Step 4: Get observation-level details
        obs_sample = adata_merged.obs.head(10).to_dict(orient="records")
        batch_counts = adata_merged.obs["batch"].value_counts().to_dict()

        # Step 5: Save the merged file
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        adata_merged.write_h5ad(
            H5AD_CONCAT,
            compression=hdf5plugin.FILTERS["zstd"],
            compression_opts=hdf5plugin.Zstd(clevel=5).filter_options,
        )

        # Response
        return {
            "merged_summary": {
                "n_cells": adata_merged.n_obs,
                "n_genes": adata_merged.n_vars,
                "batches": batch_counts,
                "output_file": H5AD_CONCAT,
                "obs_sample": obs_sample,
            }
        }
    except Exception as e:
        return {"error": str(e)}
