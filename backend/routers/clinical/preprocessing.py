from fastapi import FastAPI, HTTPException, APIRouter
from pydantic import BaseModel
import pandas as pd
import os
import logging

# Initialize FastAPI app and router

router = APIRouter()

# Constants for preprocessing
TAG_VALUES = {
    0: 0,
    1: 1,
    2: 1,
    3: 1,
    4: 1,
}

AGE_VALS_NUM = {
    "Young": 0,
    "Middle-Aged": 1,
    "Senior": 2,
}

CHOL_VALS_NUM = {
    "risk": 1,
    "normal": 0,
}

bins = [0, 30, 60, 80]
labels = ['Young', 'Middle-Aged', 'Senior']

# Configure logger
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Dataset file path
DATASET_FILE = "routers/clinical/uci_hf_df.csv"  # Path to the dataset file

# Helper function for categorical conversion
def convert_to_categorical(df, column_name, threshold):
    if column_name not in df.columns:
        raise ValueError(f"The column '{column_name}' is not in the dataframe.")
    df[column_name] = df[column_name].apply(lambda x: 'normal' if int(x) > threshold else 'risk')
    return df

# Load dataset
def load_dataset():
    if not os.path.exists(DATASET_FILE):
        logger.error(f"Dataset file '{DATASET_FILE}' not found.")
        raise FileNotFoundError(f"Dataset file '{DATASET_FILE}' not found.")
    try:
        logger.info(f"Loading dataset from file: {os.path.abspath(DATASET_FILE)}")
        return pd.read_csv(DATASET_FILE)
    except Exception as e:
        logger.error(f"Error loading dataset: {str(e)}")
        raise RuntimeError("Failed to load dataset.")

# Preprocess dataset
def preprocess_dataset(df):
    try:
        logger.info("Preprocessing dataset.")

        # Replace target values
        if 'tag' in df.columns:
            df['tag'] = df['tag'].replace(TAG_VALUES)

        # Convert "chol" column to categorical
        if 'chol' in df.columns:
            df = convert_to_categorical(df, 'chol', 260)

        # Bin "age" column and map to numerical values
        if 'age' in df.columns:
            df['age'] = pd.cut(df['age'], bins=bins, labels=labels, include_lowest=True)
            df['age'] = df['age'].replace(AGE_VALS_NUM)

        # Convert categorical "chol" to numerical
        if 'chol' in df.columns:
            df['chol'] = df['chol'].replace(CHOL_VALS_NUM)

        logger.info("Dataset preprocessing completed.")
        return df
    except Exception as e:
        logger.error(f"Error during preprocessing: {str(e)}")
        raise RuntimeError("Failed to preprocess dataset.")

# Load and preprocess dataset at startup
try:
    raw_dataframe = load_dataset()
    processed_dataframe = preprocess_dataset(raw_dataframe.copy())
except Exception as e:
    logger.error(f"Critical error during dataset loading or preprocessing: {str(e)}")
    raw_dataframe = None
    processed_dataframe = None

class ProcessedData(BaseModel):
    data: list
    columns: list

@router.get("/heart/dataset/raw", response_model=ProcessedData)
def get_raw_dataset():
    """Endpoint to fetch the raw Heart Disease dataset."""
    if raw_dataframe is None:
        raise HTTPException(status_code=500, detail="Raw dataset not available.")
    try:
        response_data = raw_dataframe.to_dict(orient="records")
        columns = list(raw_dataframe.columns)
        return ProcessedData(data=response_data, columns=columns)
    except Exception as e:
        logger.error(f"Error fetching raw dataset: {str(e)}")
        raise HTTPException(status_code=500, detail="Error fetching raw dataset.")

@router.get("/heart/dataset/processed", response_model=ProcessedData)
def get_processed_dataset():
    """Endpoint to fetch the processed Heart Disease dataset."""
    if processed_dataframe is None:
        raise HTTPException(status_code=500, detail="Processed dataset not available.")
    try:
        response_data = processed_dataframe.to_dict(orient="records")
        columns = list(processed_dataframe.columns)
        return ProcessedData(data=response_data, columns=columns)
    except Exception as e:
        logger.error(f"Error fetching processed dataset: {str(e)}")
        raise HTTPException(status_code=500, detail="Error fetching processed dataset.")

@router.get("/health")
def health_check():
    """Health check endpoint."""
    return {"status": "API is running smoothly."}

#