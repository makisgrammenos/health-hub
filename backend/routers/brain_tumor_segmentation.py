from fastapi import APIRouter, UploadFile, File, HTTPException, BackgroundTasks
from fastapi.responses import FileResponse
from models.brain_tumor_segmentation.model import BrainTumorSegmentationModel
from models.brain_tumor_segmentation.preprocessing import get_preprocessing_transforms
from models.brain_tumor_segmentation.utils import save_nifti
from models.brain_tumor_segmentation.config import MODEL_CHECKPOINT_PATH
import tempfile
import torch
import os

router = APIRouter()
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model = BrainTumorSegmentationModel(MODEL_CHECKPOINT_PATH, device)
modalities = ["T1c", "T1", "T2", "FLAIR"]

@router.post("/segment")
async def segment_mri(
    background_tasks: BackgroundTasks,

    T1c: UploadFile = File(...),
    T1: UploadFile = File(...),
    T2: UploadFile = File(...),
    FLAIR: UploadFile = File(...),
):
    try:
        # Save uploaded files to temporary paths
        temp_files = {}
        for modality, file in zip(modalities, [T1c, T1, T2, FLAIR]):
            suffix = os.path.splitext(file.filename)[1]
            with tempfile.NamedTemporaryFile(delete=False, suffix=suffix) as tmp:
                tmp.write(await file.read())
                temp_files[modality] = tmp.name

        # Preprocess images
        transforms = get_preprocessing_transforms(modalities)
        data = {modality: temp_files[modality] for modality in modalities}
        preprocessed_data = transforms(data)

        # Convert to PyTorch tensor
        input_tensor = preprocessed_data["image"].unsqueeze(0)

        # Perform segmentation
        prediction = model.predict(input_tensor)

        # Save segmentation result to NIfTI
        output_file = tempfile.NamedTemporaryFile(delete=False, suffix=".nii.gz")
        output_path = output_file.name
        save_nifti(prediction, output_path)

        # Clean up input temporary files
        for file_path in temp_files.values():
            os.unlink(file_path)

        # Schedule deletion of the output file after response is sent
        background_tasks.add_task(os.unlink, output_path)

        # Return the NIfTI file as a FileResponse
        return FileResponse(
            path=output_path,
            media_type="application/gzip",
            filename="segmentation.nii.gz"
        )

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error during segmentation: {e}")
