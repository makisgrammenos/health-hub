from fastapi import APIRouter, UploadFile, File, HTTPException
from fastapi.responses import JSONResponse
from models.skin_cancer.model import load_model, get_image_transform
from PIL import Image
import torch
import io

# Initialize router
router = APIRouter()

# Define device
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

# Load the model
model = load_model('models/skin_cancer/model.pth', device)
image_transform = get_image_transform()

@router.post("/predict")
async def predict(file: UploadFile = File(...)):
    try:
        # Load the image from the request
        image = Image.open(io.BytesIO(await file.read())).convert('RGB')

        # Preprocess the image
        transformed_image = image_transform(image).unsqueeze(0).to(device)

        # Perform inference
        with torch.no_grad():
            output = model(transformed_image)
            probability = torch.sigmoid(output).item()

            # Calculate probabilities for both classes
            probabilities = {
                "Benign": 1 - probability,
                "Malignant": probability
            }

            # Determine prediction
            prediction = "Malignant" if probabilities["Malignant"] > 0.5 else "Benign"

        # Return the result
        return JSONResponse({
            "prediction": prediction,
            "probabilities": probabilities
        })
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error processing image: {str(e)}")
