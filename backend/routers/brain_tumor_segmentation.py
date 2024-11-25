from fastapi import APIRouter, UploadFile, File, HTTPException
from models.brain_tumor_segmentation.model import load_model, segment_brain_tumor
from utils.file_handling import save_uploaded_files
import tempfile
import shutil

router = APIRouter()

# Load the model at startup
model = load_model()

@router.post("/segment")
async def segment_images(
    T1c: UploadFile = File(...),
    T1: UploadFile = File(...),
    T2: UploadFile = File(...),
    FLAIR: UploadFile = File(...),
):
    """
    Endpoint to perform brain tumor segmentation. Requires four modalities (T1c, T1, T2, FLAIR).
    """
    if not model:
        raise HTTPException(status_code=500, detail="Model not loaded. Try again later.")

    # Save uploaded files temporarily
    try:
        temp_dir = tempfile.mkdtemp()
        file_paths = save_uploaded_files(
            {"T1c": T1c, "T1": T1, "T2": T2, "FLAIR": FLAIR}, temp_dir
        )

        # Perform segmentation
        result = segment_brain_tumor(file_paths, model)

        # Clean up temporary files
        shutil.rmtree(temp_dir)

        return result

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error during segmentation: {str(e)}")
