from fastapi import FastAPI, HTTPException, APIRouter
from fastapi.responses import JSONResponse
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import umap

router = APIRouter()

tag_vals = {
    0: "control",
    1: "case",
}

# Load dataset
try:
    df = pd.read_csv("routers/clinical/uci_hf_df.csv")
    tag = pd.DataFrame(df["tag"], columns=["tag"])
    df = df.drop(["tag"], axis=1)
    tag = tag.replace({"tag": tag_vals})
except Exception as e:
    raise RuntimeError(f"Failed to load dataset: {str(e)}")

@router.get("/pca_2d")
def get_pca_2d():
    try:
        pca_2d = PCA(n_components=2)
        principal_components_2d = pca_2d.fit_transform(df)
        principal_df_2d = pd.DataFrame(data=principal_components_2d, columns=["PC 1", "PC 2"])
        principal_df_2d = pd.concat([principal_df_2d, tag], axis=1)
        return JSONResponse(content=principal_df_2d.to_dict(orient="records"))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error computing PCA 2D: {str(e)}")

@router.get("/umap_2d")
def get_umap_2d():
    try:
        umap_2d = umap.UMAP(n_components=2, random_state=42)
        umap_data_2d = umap_2d.fit_transform(df)
        umap_df_2d = pd.DataFrame(umap_data_2d, columns=["UMAP 1", "UMAP 2"])
        umap_df_2d = pd.concat([umap_df_2d, tag], axis=1)
        return JSONResponse(content=umap_df_2d.to_dict(orient="records"))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error computing UMAP 2D: {str(e)}")

@router.get("/tsne_2d")
def get_tsne_2d():
    try:
        tsne_2d = TSNE(n_components=2, random_state=42)
        tsne_data_2d = tsne_2d.fit_transform(df)
        tsne_df_2d = pd.DataFrame(tsne_data_2d, columns=["T-SNE 1", "T-SNE 2"])
        tsne_df_2d = pd.concat([tsne_df_2d, tag], axis=1)
        return JSONResponse(content=tsne_df_2d.to_dict(orient="records"))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error computing T-SNE 2D: {str(e)}")

@router.get("/pca_3d")
def get_pca_3d():
    try:
        pca_3d = PCA(n_components=3)
        principal_components_3d = pca_3d.fit_transform(df)
        principal_df_3d = pd.DataFrame(data=principal_components_3d, columns=["PC 1", "PC 2", "PC 3"])
        principal_df_3d = pd.concat([principal_df_3d, tag], axis=1)
        return JSONResponse(content=principal_df_3d.to_dict(orient="records"))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error computing PCA 3D: {str(e)}")

@router.get("/umap_3d")
def get_umap_3d():
    try:
        umap_3d = umap.UMAP(n_components=3, random_state=42)
        umap_result_3d = umap_3d.fit_transform(df)
        umap_df_3d = pd.DataFrame(umap_result_3d, columns=["UMAP 1", "UMAP 2", "UMAP 3"])
        umap_df_3d = pd.concat([umap_df_3d, tag], axis=1)
        return JSONResponse(content=umap_df_3d.to_dict(orient="records"))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error computing UMAP 3D: {str(e)}")

@router.get("/tsne_3d")
def get_tsne_3d():
    try:
        tsne_3d = TSNE(n_components=3, random_state=42)
        tsne_data_3d = tsne_3d.fit_transform(df)
        tsne_df_3d = pd.DataFrame(tsne_data_3d, columns=["T-SNE 1", "T-SNE 2", "T-SNE 3"])
        tsne_df_3d = pd.concat([tsne_df_3d, tag], axis=1)
        return JSONResponse(content=tsne_df_3d.to_dict(orient="records"))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error computing T-SNE 3D: {str(e)}")


