from fastapi import APIRouter, UploadFile, File, HTTPException
from PIL import Image
from io import BytesIO
from models.covidnext50.loader import COVIDNext50Model
import os

router = APIRouter()

# Initialize the model
model_path =  "models/covidnext50/checkpoint/model.pth"
print(os.getcwd()) 
model = COVIDNext50Model(model_path)

@router.post("/predict-covid")
async def predict_covid(file: UploadFile = File(...)):
    """
    Predict COVID presence in a chest X-ray image.
    """
    try:
        # Read image
        contents = await file.read()
        image = Image.open(BytesIO(contents)).convert("RGB")

        # Predict
        prediction,probabilities = model.predict(image)

        return {"message": "Prediction successful", "prediction": prediction, "probabilities": probabilities}
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Prediction failed: {e}")
