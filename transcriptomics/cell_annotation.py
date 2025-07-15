from fastapi import APIRouter, HTTPException
from fastapi.responses import JSONResponse
import scanpy as sc
import decoupler as dc
import pandas as pd
import numpy as np
import hdf5plugin
import matplotlib.pyplot as plt
import seaborn as sns
from io import BytesIO
import base64
import os

router = APIRouter()

# Paths
H5AD_PATH = "routers/transcriptomics/data/4.H5AD_integrated/GSE161872_scanorama_v1.h5ad"
MARKER_GENES_PATH = "routers/transcriptomics/data/marker_genes/mouse_adipose_tissue_cellmarker2_act_combination.csv"
OUTPUT_H5AD_PATH = "routers/transcriptomics/data/5.H5AD_annotated_clusters/GSE161872_scanorama_integration_decoupleR_anno_v1.h5ad"

CLUSTER_RESOLUTION = 0.25


def create_marker_dict(df_cellmarker, cell_name="cell_name", marker="Symbol"):
    cell_markers = {}
    for i in set(df_cellmarker[cell_name].values.tolist()):
        genes = df_cellmarker[df_cellmarker[cell_name] == i][marker].values.tolist()
        cell_markers[i] = genes
    return cell_markers


def create_melted_df(score_df, ctype_lst):
    dfs = [score_df[score_df['cluster'] == str(i)].assign(dataset=f'dataset_{i}') for i in range(len(score_df["cluster"].value_counts()))]
    final_df = pd.concat(dfs, ignore_index=True)
    final_df_melted = final_df.melt(id_vars=['dataset', 'cluster'], 
                                    value_vars=ctype_lst,
                                    var_name='cell_type',
                                    value_name='score')
    return final_df_melted


@router.post("/ora_analysis")
async def perform_ora_analysis():
    try:
        # Load the data
        if not os.path.exists(H5AD_PATH):
            raise HTTPException(status_code=404, detail=f"{H5AD_PATH} does not exist.")

        adata = sc.read_h5ad(H5AD_PATH)

        # Handle NaN values in the data matrix
        if np.isnan(adata.X).any():
            raise HTTPException(
                status_code=422,
                detail="Input data contains NaN values. Please ensure preprocessing removes or imputes missing values."
            )

        # Perform Leiden clustering
        sc.tl.leiden(adata, resolution=CLUSTER_RESOLUTION, key_added="leiden")

        # Read marker genes
        cell_markers = pd.read_csv(MARKER_GENES_PATH).drop_duplicates()

        # Run ORA
        dc.run_ora(
            mat=adata,
            net=cell_markers,
            source="cell_name",
            target="Symbol",
            min_n=3,
            verbose=True,
            use_raw=False,
        )

        # Get ORA scores
        acts = dc.get_acts(adata, obsm_key="ora_estimate")
        acts_v = acts.X.ravel()
        max_e = np.nanmax(acts_v[np.isfinite(acts_v)])
        acts.X[~np.isfinite(acts.X)] = max_e

        # Create cell-type list and ORA-score dataframe
        score_df = adata.obsm['ora_estimate']
        ctype_lst = list(score_df.columns)
        score_df["cluster"] = adata.obs["leiden"]

        # Melted dataframe for plotting
        melted_df = create_melted_df(score_df, ctype_lst)

        # Plot Leiden clustering
        plt.figure()
        sc.pl.umap(adata, color='leiden', add_outline=True, legend_loc="on data",
                   legend_fontsize=12, legend_fontoutline=2, frameon=False, title='Leiden Clustering', show=False)
        leiden_plot = BytesIO()
        plt.savefig(leiden_plot, format="png", dpi=300, bbox_inches="tight")
        plt.close()

        # Generate violin plot
        plt.figure(figsize=(16, 6))
        sns.violinplot(x='cluster', y='score', hue='cell_type', data=melted_df, 
                       inner="quartile", palette="muted", scale="count")
        plt.title("ORA Scores Across Clusters")
        plt.ylabel("ORA Score")
        plt.xlabel("Leiden Cluster")
        plt.legend(title="Cell Type", bbox_to_anchor=(1.01, 1), loc="upper left")
        violin_plot = BytesIO()
        plt.savefig(violin_plot, format="png", dpi=300, bbox_inches="tight")
        plt.close()

        # Save updated dataset
        # adata.write_h5ad(OUTPUT_H5AD_PATH, compression=hdf5plugin.FILTERS["zstd"])

        return JSONResponse(
            content={
                "message": "ORA analysis completed successfully.",
                "output_file": OUTPUT_H5AD_PATH,
                "leiden_plot": base64.b64encode(leiden_plot.getvalue()).decode("utf-8"),
                "violin_plot": base64.b64encode(violin_plot.getvalue()).decode("utf-8"),
            }
        )

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"ORA analysis failed: {str(e)}")
