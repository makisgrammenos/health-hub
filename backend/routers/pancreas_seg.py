# models/pancreas_seg/router.py
import os
import uuid
import shutil
import tempfile
import logging

from fastapi import APIRouter, UploadFile, File, HTTPException, BackgroundTasks
from fastapi.responses import FileResponse

from models.pancreas_seg.model import predict_single_nifti

logger = logging.getLogger(__name__)
router = APIRouter()

def _is_nifti(name: str) -> bool:
    n = name.lower()
    return n.endswith(".nii") or n.endswith(".nii.gz")

@router.post("/segment")
async def segment_ct(background_tasks: BackgroundTasks, file: UploadFile = File(...)):
    if not _is_nifti(file.filename):
        raise HTTPException(400, "Upload a NIfTI volume (.nii or .nii.gz).")

    # persist upload
    suffix = ".nii.gz" if file.filename.lower().endswith(".nii.gz") else ".nii"
    with tempfile.NamedTemporaryFile(delete=False, suffix=suffix) as tmp:
        shutil.copyfileobj(file.file, tmp)
        tmp_in = tmp.name

    try:
        mask_path, overlays, volumes = await predict_single_nifti(tmp_in)
        # clean the upload once we’re done
        background_tasks.add_task(os.unlink, tmp_in)
        # if you want the mask to remain downloadable for a while, DO NOT delete it here.
        # otherwise, uncomment the line below to delete after response is sent:
        # background_tasks.add_task(os.unlink, mask_path)

        filename = os.path.basename(mask_path)
        return {
            "segmented_slices": overlays,                        # base64 PNGs
            "download_url": f"/imaging/pancreas/download/{filename}",
            "volumes": volumes,
            "status": "success",
        }
    except Exception as e:
        logger.exception("Pancreas segmentation failed")
        try:
            os.unlink(tmp_in)
        except Exception:
            pass
        raise HTTPException(500, f"Error during segmentation: {e}") from e

@router.get("/download/{filename}")
async def download_segmentation(filename: str):
    path = os.path.join(tempfile.gettempdir(), filename)
    if not os.path.exists(path):
        raise HTTPException(404, "File not found or expired.")
    return FileResponse(path, media_type="application/gzip", filename="pancreas_segmentation.nii.gz")
