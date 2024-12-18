from fastapi import APIRouter, UploadFile, File, HTTPException
from models.breast_cancer.model import BreastCancerModel
from typing import List
from pydantic import BaseModel
from PIL import Image
import zipfile
import os
import tempfile

router = APIRouter()

# Initialize the model
model = BreastCancerModel("models/breast_cancer/checkpoint/model.pth")

# Response models
class PredictionResponse(BaseModel):
    filename: str
    label: str
    probabilities: dict

class BatchPredictionResponse(BaseModel):
    message: str
    results: List[PredictionResponse]

@router.post("/predict")
async def predict_single_image(file: UploadFile = File(...)):
    try:
        image = Image.open(file.file).convert("RGB")
        label, probabilities = model.predict(image)
        print(label, probabilities)
        if label == 0:
            label = "No Invasive Ductal Carcinoma"
        elif label == 1:
            label = "Invasive Ductal Carcinoma"
        probabilities = {
            "No Invasive Ductal Carcinoma": probabilities[0][0],
            "Invasive Ductal Carcinoma": probabilities[0][1]
        }
        
        return {"filename": file.filename, "result": label, "probabilities": probabilities}
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Prediction failed: {e}")

@router.post("/predict-zip")
async def predict_zip_file(file: UploadFile = File(...)):
    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            zip_path = os.path.join(tmpdir, "uploaded.zip")
            with open(zip_path, "wb") as temp_zip:
                temp_zip.write(file.file.read())

            with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                zip_ref.extractall(tmpdir)

            results = []
            for root, _, files in os.walk(tmpdir):
                for img_file in files:
                    try:
                        img_path = os.path.join(root, img_file)
                        image = Image.open(img_path).convert("RGB")
                        label, probabilities = model.predict(image)
                        response ={
                            "filename": img_file,
                            "result": label,
                            "probabilities": probabilities
                        }
                        results.append(response)
                    except Exception as e:
                        results.append({"filename": img_file, "error": str(e)})

        return  {"message": "Batch processing completed", "results": results}
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Batch processing failed: {e}")
