from fastapi import APIRouter, HTTPException
import scanpy as sc
import numpy as np
import os
import hdf5plugin
import scanorama
import warnings

warnings.filterwarnings("ignore")

router = APIRouter()

# Paths
H5AD = "routers/transcriptomics/data/3.H5AD_concatenated/GSE161872_v1.h5ad"
H5AD_INTEGRATED = "routers/transcriptomics/data/4.H5AD_integrated/GSE161872_scanorama_v1.h5ad"

def data_preprocessor(adata):
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    sc.pp.scale(adata, max_value=10)
    sc.pp.pca(adata)
    return adata

@router.post("/integrate")
async def integrate_data():
    try:
        # Read and preprocess
        merged_adata = sc.read_h5ad(H5AD)
        merged_adata = data_preprocessor(merged_adata)

        # Split batches
        batches = merged_adata.obs['batch'].cat.categories.tolist()
        alldata = {batch: merged_adata[merged_adata.obs['batch'] == batch,] for batch in batches}

        # Integrate with Scanorama
        adatas = list(alldata.values())
        scanorama.integrate_scanpy(adatas, dimred=50)
        scanorama_int = [ad.obsm['X_scanorama'] for ad in adatas]
        all_s = np.concatenate(scanorama_int)
        adata_all_scanorama = merged_adata.copy()
        adata_all_scanorama.obsm["Scanorama"] = all_s

        # Calculate UMAP
        sc.pp.neighbors(adata_all_scanorama, use_rep="Scanorama")
        sc.tl.umap(adata_all_scanorama)

        # Save integrated file
        adata_all_scanorama.write_h5ad(
            H5AD_INTEGRATED,
            compression=hdf5plugin.FILTERS["zstd"]
        )

        return {
            "message": "Integration completed successfully.",
            "output_file": H5AD_INTEGRATED,
            "batches": {batch: adata_all_scanorama[adata_all_scanorama.obs['batch'] == batch].shape[0] for batch in batches}
        }

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Integration failed: {str(e)}")
