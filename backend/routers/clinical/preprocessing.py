from fastapi import File, UploadFile, HTTPException, Form, APIRouter
from fastapi.responses import JSONResponse
from fastapi.middleware.cors import CORSMiddleware
import pandas as pd
import io

router = APIRouter()

def convert_to_categorical(df, column_name, threshold):
    if column_name not in df.columns:
        raise ValueError(f"The column '{column_name}' is not in the dataframe.")

    df[column_name] = df[column_name].apply(lambda x: 'normal' if int(x) > threshold else 'risk')
    return df

def bin_column(df, column_name, bins, labels):
    if column_name not in df.columns:
        raise ValueError(f"The column '{column_name}' is not in the dataframe.")

    df[column_name] = pd.cut(df[column_name], bins=bins, labels=labels, include_lowest=True)
    return df

def fill_missing_values(df, column_name, value):
    if column_name not in df.columns:
        raise ValueError(f"The column '{column_name}' is not in the dataframe.")

    df[column_name] = df[column_name].fillna(value)
    return df

def preprocess_data(df, options):
    # Handle missing values
    if options.get("drop_missing", False):
        df = df.dropna()
    elif options.get("fill_missing"):
        fill_options = options.get("fill_missing")
        for column, value in fill_options.items():
            df = fill_missing_values(df, column, value)

    # Convert to categorical if requested
    if options.get("categorical_column") and options.get("threshold") is not None:
        column = options["categorical_column"]
        threshold = options["threshold"]
        df = convert_to_categorical(df, column, threshold)
        df[column] = df[column].replace({"risk": 1, "normal": 0})

    # Binning if requested
    if options.get("bin_column") and options.get("bins") and options.get("labels"):
        column = options["bin_column"]
        bins = options["bins"]
        labels = options["labels"]
        df = bin_column(df, column, bins, labels)
        df[column] = df[column].replace({"Young": 0, "Middle-Aged": 1, "Senior": 2})

    return df

@router.post("/basic-processing/upload")
async def upload_dataset(
    file: UploadFile = File(...),
    drop_missing: bool = Form(False),
    fill_column: str = Form(None),
    fill_value: str = Form(None),
    categorical_column: str = Form(None),
    threshold: float = Form(None),
    bin_column: str = Form(None),
    bins: str = Form(None),  # JSON string representation of list
    labels: str = Form(None)  # JSON string representation of list
):
    try:
        # Read uploaded file
        contents = await file.read()
        df = pd.read_csv(io.StringIO(contents.decode("utf-8")))

        # Parse processing options
        options = {
            "drop_missing": drop_missing,
            "fill_missing": {fill_column: fill_value} if fill_column and fill_value else None,
            "categorical_column": categorical_column,
            "threshold": threshold,
            "bin_column": bin_column,
            "bins": eval(bins) if bins else None,
            "labels": eval(labels) if labels else None,
        }

        # Preprocess dataset
        processed_df = preprocess_data(df, options)

        # Convert the processed dataframe to JSON
        processed_json = processed_df.to_dict(orient="records")
        return JSONResponse(content={"message": "Dataset processed successfully!", "data": processed_json})

    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

@router.post("/preview-columns")
async def preview_columns(file: UploadFile = File(...)):
    try:
        # Read uploaded file
        contents = await file.read()
        df = pd.read_csv(io.StringIO(contents.decode("utf-8")))

        # Get column names
        columns = df.columns.tolist()
        return JSONResponse(content={"columns": columns})

    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

@router.get("/example/processed-dataset")
def get_example_dataset():
    # Provide an example preprocessed dataset for testing
    example_data = {
        "age": [0, 1, 2],
        "chol": [1, 0, 1],
        "other_column": [10, 20, 30],
    }
    return JSONResponse(content=example_data)

# Run the application
# uvicorn script_name:app --reload
