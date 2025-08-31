from fastapi import APIRouter, UploadFile, File, HTTPException
from pydantic import BaseModel
from typing import Dict
from PIL import Image
import io
import os

from models.breast_density.model import BiradsDensityModel

router = APIRouter()

# ENV toggles (optional)
CHECKPOINT_PATH = os.getenv("BIRADS_CHECKPOINT", "models/breast_density/model.pt")
DEBUG_FLAG = bool(int(os.getenv("BIRADS_DEBUG", "0")))

# Initialize once (same weights & preprocessing as CLI script)
MODEL = BiradsDensityModel(checkpoint_path=CHECKPOINT_PATH, debug=DEBUG_FLAG)


class PredictionResponse(BaseModel):
    filename: str
    prediction_index: int
    prediction_label: str
    probabilities: Dict[str, float]


@router.post("/predict", response_model=PredictionResponse)
async def predict_single_image(file: UploadFile = File(...)):
    # Same simple checks you already used
    if file.content_type not in {"image/jpeg", "image/png", "image/tiff", "image/bmp"}:
        raise HTTPException(status_code=400, detail=f"Unsupported content type: {file.content_type}")

    raw = await file.read()
    if not raw:
        raise HTTPException(status_code=400, detail="Empty file.")

    # Decode with PIL
    try:
        pil = Image.open(io.BytesIO(raw))
    except Exception as e:
        raise HTTPException(status_code=400, detail=f"Invalid image: {e}")

    try:
        pred_idx, probs = MODEL.predict(pil)
        return MODEL.format_response(file.filename, pred_idx, probs)
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Prediction failed: {e}")
