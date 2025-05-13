from fastapi import APIRouter, HTTPException
from fastapi.responses import JSONResponse
import pandas as pd
from pycaret.classification import setup, compare_models, predict_model, pull

router = APIRouter()

@router.get("/classification")
def run_pycaret_classification():
    try:
        # Load the dataset
        df = pd.read_csv("routers/clinical/uci_hf_df.csv")

        # Setup the PyCaret environment
        s = setup(df, target='tag', session_id=123,  verbose=False)

        # Compare models and select the best
        best_model = compare_models()

        # Fetch performance metrics
        metrics = pull()

        # Generate predictions on the holdout dataset
        holdout_predictions = predict_model(best_model)

        # Return metrics and predictions
        return JSONResponse(content={
            "metrics": metrics.to_dict(orient="records"),
            "predictions": holdout_predictions.to_dict(orient="records")
        })
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error running PyCaret classification: {str(e)}")
