from fastapi import FastAPI, HTTPException, APIRouter
from fastapi.responses import JSONResponse
import pandas as pd
from sklearn.preprocessing import StandardScaler

router = APIRouter()

TAG_VALS_CAT = {
    1: "case",
    0: "control",
}

TAG_VALS_NUM = {
    "case": 1,
    "control": 0,
}

# Load dataset
try:
    df = pd.read_csv("routers/clinical/uci_hf_df.csv")
    df["tag"] = df["tag"].replace(TAG_VALS_CAT)
except Exception as e:
    raise RuntimeError(f"Failed to load dataset: {str(e)}")

# Utility function to sanitize data
def sanitize_dataframe(df):
    """Replace non-JSON-compliant values in the DataFrame."""
    return df.replace([float('inf'), float('-inf')], None).fillna(0)

# Sanitize the dataset globally
df = sanitize_dataframe(df)

@router.get("/category_counts")
def get_category_counts():
    try:
        # Count occurrences of each category
        category_counts = df["tag"].value_counts().reindex(["control", "case"]).fillna(0).to_dict()
        return JSONResponse(content={"category_counts": category_counts})
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error computing category counts: {str(e)}")

@router.get("/standardized_data")
def get_standardized_data():
    try:
        # Standardize numerical data
        tag_df = df[["tag"]].copy()
        df_numerical = df.drop("tag", axis=1)
        scaler = StandardScaler()
        df_scaled_numerical = pd.DataFrame(
            scaler.fit_transform(df_numerical),
            columns=df_numerical.columns,
            index=df_numerical.index
        )
        df_scaled = pd.concat([df_scaled_numerical, tag_df], axis=1)

        return JSONResponse(content=df_scaled.to_dict(orient="records"))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error standardizing data: {str(e)}")

@router.get("/melted_data")
def get_melted_data():
    try:
        # Melt scaled data
        tag_df = df[["tag"]].copy()
        df_numerical = df.drop("tag", axis=1)
        scaler = StandardScaler()
        df_scaled_numerical = pd.DataFrame(
            scaler.fit_transform(df_numerical),
            columns=df_numerical.columns,
            index=df_numerical.index
        )
        df_scaled = pd.concat([df_scaled_numerical, tag_df], axis=1)
        melted_df = pd.melt(
            df_scaled,
            id_vars="tag",
            value_vars=[col for col in df_scaled.columns if col != "condition"],
            var_name="feature",
            value_name="value"
        )

        return JSONResponse(content=melted_df.to_dict(orient="records"))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error melting data: {str(e)}")

@router.get("/correlation_matrix")
def get_correlation_matrix():
    try:
        df["tag"] = df["tag"].replace(TAG_VALS_NUM)
        corr_matrix = df.corr(method="pearson")
        return JSONResponse(content=corr_matrix.to_dict())
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error computing correlation matrix: {str(e)}")

app = FastAPI()
app.include_router(router, prefix="/demo", tags=["Data Visualization Demo"])
