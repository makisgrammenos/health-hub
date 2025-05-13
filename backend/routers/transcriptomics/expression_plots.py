import base64
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
import os
import scanpy as sc
import matplotlib.pyplot as plt

router = APIRouter()

# Paths
H5AD_PATH = "routers/transcriptomics/data/5.H5AD_annotated_clusters/GSE161872_scanorama_integration_decoupleR_anno_v1.h5ad"
FIGURE_PATH = "./figures/umap_expression_plot.png"

GENES = ['Man2b1', 'Pld3', 'Evi2a', 'Ar', 'Col6a1', 'Gstm2']

# Utility function
def plot_gene_expr(adata, gene_lst, rows, columns):
    fig, axs = plt.subplots(rows, columns, figsize=(7, 4))
    axs = axs.flatten()

    for i, gene in enumerate(gene_lst):
        if i < len(axs):
            sc.pl.umap(adata, color=gene, add_outline=True, legend_loc='on data',
                       legend_fontsize=12, legend_fontoutline=2, frameon=True,
                       title=f'{gene}', palette='Set1', ax=axs[i], show=False)

    os.makedirs(os.path.dirname(FIGURE_PATH), exist_ok=True)
    plt.tight_layout()
    plt.savefig(FIGURE_PATH, dpi=600, bbox_inches='tight')
    plt.close(fig)

# Endpoint
@router.get("/umap-expression")
async def get_umap_expression():
    try:
        if not os.path.exists(H5AD_PATH):
            raise HTTPException(status_code=404, detail="H5AD file not found.")

        adata = sc.read_h5ad(H5AD_PATH)
        plot_gene_expr(adata, GENES, 2, 3)

        if not os.path.exists(FIGURE_PATH):
            raise HTTPException(status_code=500, detail="Failed to generate the UMAP plot.")

        # Read and encode the image in Base64
        with open(FIGURE_PATH, "rb") as img_file:
            encoded_string = base64.b64encode(img_file.read()).decode("utf-8")

        return {"image_base64": encoded_string}

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error: {str(e)}")
