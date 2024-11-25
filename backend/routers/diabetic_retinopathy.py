from fastapi import APIRouter, UploadFile, File, HTTPException
from models.diabetic_retinopathy.model import predict_single, predict_batch
from utils.file_handling import extract_zip

router = APIRouter()

@router.post("/predict")
async def predict_image(file: UploadFile = File(...)):
    if file.content_type not in ["image/jpeg", "image/png"]:
        raise HTTPException(status_code=400, detail="Invalid file type. Please upload a JPEG or PNG image.")
    try:
        prediction = await predict_single(file)
        return prediction
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error during prediction: {str(e)}")

@router.post("/batch-predict")
async def batch_predict(file: UploadFile = File(...)):
    if file.content_type != "application/zip":
        raise HTTPException(status_code=400, detail="Invalid file type. Please upload a ZIP file.")
    try:
        temp_dir = extract_zip(file)
        predictions = await predict_batch(temp_dir)
        return {"predictions": predictions}
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error during batch prediction: {str(e)}")
