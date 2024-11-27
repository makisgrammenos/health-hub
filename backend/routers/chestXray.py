from fastapi import APIRouter, UploadFile, File, HTTPException, BackgroundTasks
from fastapi.responses import JSONResponse
from models.chest_xray.model import ChestXRayModel
from models.chest_xray.preprocessing import preprocess_image
import tempfile
import zipfile
import os
import numpy as np
from PIL import Image

router = APIRouter()
model = ChestXRayModel()  # Load model instance


@router.post("/predict")
async def classify_single_image(file: UploadFile = File(...)):
    """
    Classify a single chest X-ray image.
    """
    try:
        # Save the uploaded file to a temporary location
        with tempfile.NamedTemporaryFile(delete=False, suffix=".tmp") as tmp:
            tmp.write(await file.read())
            temp_filepath = tmp.name

        # Read the image
        if file.content_type == "application/dicom" or file.filename.lower().endswith(".dcm"):
            from pydicom import dcmread
            dicom_data = dcmread(temp_filepath)
            img = dicom_data.pixel_array
            image_display = Image.fromarray(img)
        else:
            img = Image.open(temp_filepath).convert("RGB")
            image_display = img
            img = np.array(img)

        # Preprocess and run inference
        img_tensor = preprocess_image(img, model.device)
        outputs = model.predict(img_tensor)
        outputs = outputs.cpu().numpy()[0]
        pathologies = model.get_pathologies()

        # Organize predictions
        predictions = [
            {"Pathology": pathology, "Probability": f"{prob * 100:.2f}%"}
            for pathology, prob in zip(pathologies, outputs)
        ]
        predictions = sorted(predictions, key=lambda x: x["Probability"], reverse=True)

        # Clean up temporary file
        os.unlink(temp_filepath)

        return JSONResponse(content={"predictions": predictions})
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error processing the image: {e}")


@router.post("/batch-prediction")
async def classify_batch_images(file: UploadFile = File(...)):
    
    """
    Classify a batch of chest X-ray images uploaded as a ZIP file.
    """
    try:
        with tempfile.TemporaryDirectory() as tempdir:
            # Save uploaded ZIP to a temporary file
            zip_path = os.path.join(tempdir, "batch.zip")
            with open(zip_path, "wb") as f:
                f.write(await file.read())

            # Extract ZIP file
            with zipfile.ZipFile(zip_path, "r") as zip_ref:
                zip_ref.extractall(tempdir)

            # Find valid image files
            valid_extensions = [".jpg", ".jpeg", ".png", ".dcm"]
            image_files = [
                os.path.join(root, file)
                for root, _, files in os.walk(tempdir)
                for file in files
                if os.path.splitext(file)[1].lower() in valid_extensions
            ]

            if not image_files:
                raise HTTPException(status_code=400, detail="No valid image files found in the ZIP.")

            results = []
            for image_path in image_files:
                try:
                    # Read and preprocess the image
                    if image_path.endswith(".dcm"):
                        from pydicom import dcmread
                        dicom_data = dcmread(image_path)
                        img = dicom_data.pixel_array
                    else:
                        img = Image.open(image_path).convert("RGB")
                        img = np.array(img)

                    img_tensor = preprocess_image(img, model.device)
                    outputs = model.predict(img_tensor)
                    outputs = outputs.cpu().numpy()[0]
                    pathologies = model.get_pathologies()

                    # Get top prediction
                    top_idx = np.argmax(outputs)
                    results.append({
                        "Image": os.path.basename(image_path),
                        "Top Pathology": pathologies[top_idx],
                        "Probability": f"{outputs[top_idx] * 100:.2f}%"
                    })
                except Exception as e:
                    results.append({"Image": os.path.basename(image_path), "Error": str(e)})

            return JSONResponse(content={"predictions": results})
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error processing the ZIP file: {e}")
