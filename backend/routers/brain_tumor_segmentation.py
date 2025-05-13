from fastapi import APIRouter, UploadFile, File, HTTPException, BackgroundTasks
from fastapi.responses import FileResponse
from models.brain_tumor_segmentation.model import BrainTumorSegmentationModel
from models.brain_tumor_segmentation.preprocessing import get_preprocessing_transforms
from models.brain_tumor_segmentation.utils import save_nifti
from models.brain_tumor_segmentation.config import MODEL_CHECKPOINT_PATH

import tempfile
import torch
import os
import nibabel as nib
import numpy as np
import base64
from PIL import Image
import io
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,ListedColormap,BoundaryNorm
import logging

# --------------------------------------------------------------------------- #
# logging & global initialisation
# --------------------------------------------------------------------------- #
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

router = APIRouter()
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

try:
    model = BrainTumorSegmentationModel(MODEL_CHECKPOINT_PATH, device)
    logger.info(f"Model loaded successfully on {device}")
except Exception as e:
    logger.error(f"Failed to load model: {e}")
    model = None

modalities = ["T1c", "T1", "T2", "FLAIR"]

# # # colourmap for overlay
# tumor_colors = [
#     (0, 0, 0, 0),      # 0 background
#     (1, 0, 0, 0.7),    # 1 necrotic / non-enhancing
#     (0, 0, 1, 0.7),    # 2 edema
#     (0, 0, 0, 0),      # 3 unused
#     (0, 1, 0, 0.7)     # 4 enhancing
# ]
# tumor_cmap = LinearSegmentedColormap.from_list("tumor_cmap", tumor_colors, N=5)
tumor_colors = [
    (0, 0, 0, 0),      # 0 background
    (0, 1, 0, 0.7),    # 1 green  → necrotic/non-enhancing  (optional)
    (0, 0, 1, 0.7),    # 2 blue   → oedema
    (1, 0, 0, 0.7),    # 3 red    → enhancing tumour **(what you called "tumor")**
]
tumor_cmap = ListedColormap(tumor_colors, name="tumor_cmap")
tumor_norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], ncolors=len(tumor_colors))

# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def create_overlay_image(mri_slice: np.ndarray, seg_slice: np.ndarray) -> str:
    """
    Return a base-64 PNG with the segmentation overlaid on the MRI slice.
    Both inputs must be 2-D arrays with identical shape.
    """
    mri_norm = (
        (mri_slice - mri_slice.min()) / (mri_slice.ptp() + 1e-8)
        if mri_slice.max() > mri_slice.min()
        else np.zeros_like(mri_slice)
    )

    plt.figure(figsize=(10, 10), dpi=100)
    plt.imshow(mri_norm, cmap="gray")
    # plt.imshow(seg_slice, cmap=tumor_cmap, alpha=0.7)
    plt.imshow(seg_slice, cmap=tumor_cmap, norm=tumor_norm, alpha=0.7)

    plt.axis("off")
    plt.tight_layout(pad=0)

    buf = io.BytesIO()
    plt.savefig(buf, format="png", bbox_inches="tight", pad_inches=0)
    plt.close()
    buf.seek(0)

    return "data:image/png;base64," + base64.b64encode(buf.read()).decode("utf-8")


# --------------------------------------------------------------------------- #
# routes
# --------------------------------------------------------------------------- #
@router.post("/segment")
async def segment_mri(
    background_tasks: BackgroundTasks,
    T1c: UploadFile = File(...),
    T1: UploadFile = File(...),
    T2: UploadFile = File(...),
    FLAIR: UploadFile = File(...),
):
    if model is None:
        raise HTTPException(
            status_code=503,
            detail="Model is not available. Please try again later.",
        )

    temp_files, output_path = {}, None
    try:
        # ------------------------------------------------------------------ #
        # 1. save incoming files to disk so nibabel can open them
        # ------------------------------------------------------------------ #
        for modality, f in zip(modalities, [T1c, T1, T2, FLAIR]):
            if not f.filename.endswith((".nii", ".nii.gz")):
                raise HTTPException(
                    status_code=400,
                    detail=f"Invalid file format for {modality}. Only .nii or .nii.gz files are supported.",
                )

            suffix = ".nii.gz" if f.filename.endswith(".nii.gz") else ".nii"
            with tempfile.NamedTemporaryFile(delete=False, suffix=suffix) as tmp:
                tmp.write(await f.read())
                temp_files[modality] = tmp.name
            logger.info(f"Saved {modality} file to {temp_files[modality]}")

        #------------------------------------------------------------------  #
        # 2. preprocessing & forward pass
        # ------------------------------------------------------------------ #
        transforms = get_preprocessing_transforms(modalities)
        data_dict = {m: temp_files[m] for m in modalities}
        preprocessed = transforms(data_dict)  # -> dict with "image"

        input_tensor = preprocessed["image"].unsqueeze(0).to(device)  # add batch dim
        pred_np = model.predict(input_tensor)  # torch Tensor

        # ------------------------------------------------------------------ #
        # 3. convert network output to label map (H,W,D)
        # ------------------------------------------------------------------ #
        # pred_np = prediction.squeeze().detach().cpu().numpy()  # <-- CONVERT HERE

        # if pred_np.ndim == 4:  # multi-class logits
        #     label_map = np.argmax(pred_np, axis=0).astype(np.uint8)
        # elif pred_np.ndim == 3:  # already single-chan
        #     label_map = pred_np.astype(np.uint8)
        # else:
        #     raise RuntimeError(f"Unexpected prediction shape: {pred_np.shape}")
        # unique, counts = np.unique(label_map, return_counts=True)
        # print(dict(zip(unique, counts)))

        wt = pred_np[0] > 0.5  # whole tumour          → label 2 (blue)
        tc = pred_np[1] > 0.5  # tumour core           → label 1 (red)
        et = pred_np[2] > 0.5  # enhancing tumour      → label 3 (green)

        label_map = np.zeros(pred_np.shape[1:], dtype=np.uint8)
        label_map[wt] = 2
        label_map[tc] = 1
        label_map[et] = 3

        # ------------------------------------------------------------------ #
        # 4. save nifti
        # ------------------------------------------------------------------ #
        reference_img = nib.load(temp_files["T1c"])
        output_file = tempfile.NamedTemporaryFile(delete=False, suffix=".nii.gz")
        output_path = output_file.name
        output_file.close()

        save_nifti(label_map, output_path, affine=reference_img.affine)  # <-- PASS NumPy

        # # ------------------------------------------------------------------ #
        # # 5. slice selection & visualisation (axial ⇒ third axis)
        # # ------------------------------------------------------------------ #
        # nz_slices = np.where(np.sum(label_map > 0, axis=(0, 1)) > 10)[0]
        # logger.info(f"Found {len(nz_slices)} slices with segmentation")
        #
        # if len(nz_slices) == 0:
        #     mid = label_map.shape[2] // 2
        #     slice_indices = [i for i in (mid - 5, mid, mid + 5) if 0 <= i < label_map.shape[2]]
        #     logger.warning("No tumour voxels detected – using middle slices")
        # else:
        #     # first, middle, last tumour slice – or all of them if < 3
        #     if len(nz_slices) >= 3:
        #         slice_indices = [nz_slices[0], nz_slices[len(nz_slices) // 2], nz_slices[-1]]
        #     else:
        #         slice_indices = nz_slices.tolist()
        # logger.info(f"Selected slices: {slice_indices}")
        #
        # reference_vol = reference_img.get_fdata()           # same shape as label_map
        # segmented_slices = []
        # for idx in slice_indices:
        #     mri_slice = reference_vol[:, :, idx]
        #     seg_slice = label_map[:, :, idx]
        #     segmented_slices.append(
        #         {"slice_index": int(idx), "image": create_overlay_image(mri_slice, seg_slice)}
        #     )
        # ------------------------------------------------------------------ #
        # 5. slice visualisation – **return ALL slices**
        # ------------------------------------------------------------------ #
        reference_vol = reference_img.get_fdata()  # same shape as label_map
        segmented_slices = []

        # iterate over every axial slice
        for idx in range(label_map.shape[2]):
            mri_slice = reference_vol[:, :, idx]
            seg_slice = label_map[:, :, idx]
            segmented_slices.append(
                {"slice_index": int(idx),
                 "image": create_overlay_image(mri_slice, seg_slice)}
            )
        logger.info(f"Prepared {len(segmented_slices)} slices for front-end")

        # ------------------------------------------------------------------ #
        # 6. statistics
        # ------------------------------------------------------------------ #
        # tumor_volumes = {
        #     "necrotic": float(np.sum(label_map == 1)),
        #     "edema": float(np.sum(label_map == 2)),
        #     "enhancing": float(np.sum(label_map == 4)),
        #     "total": float(np.sum(label_map > 0)),
        # }
        tumor_volumes = {
            "necrotic": float(np.sum(label_map == 1)),
            "edema": float(np.sum(label_map == 2)),
            "enhancing": float(np.sum(label_map == 3)),
            "total": float(np.sum(label_map > 0)),
        }

        # ------------------------------------------------------------------ #
        # 7. schedule temp-file cleanup and respond
        # ------------------------------------------------------------------ #
        for p in temp_files.values():
            background_tasks.add_task(os.unlink, p)
        background_tasks.add_task(os.unlink, output_path)

        return {
            "segmented_slices": segmented_slices,
            "download_url": f"/imaging/brain-tumor/download/{os.path.basename(output_path)}",
            "tumor_volumes": tumor_volumes,
            "status": "success",
        }

    except Exception as exc:
        logger.error("Error during segmentation", exc_info=True)

        # tidy up anything we created
        for p in temp_files.values():
            if os.path.exists(p):
                try:
                    os.unlink(p)
                except Exception:
                    pass
        if output_path and os.path.exists(output_path):
            try:
                os.unlink(output_path)
            except Exception:
                pass

        raise HTTPException(status_code=500, detail=f"Error during segmentation: {exc}") from exc


@router.get("/download/{filename}")
async def download_segmentation(filename: str):
    path = os.path.join(tempfile.gettempdir(), filename)
    if not os.path.exists(path):
        raise HTTPException(status_code=404, detail="File not found or has expired")

    try:
        return FileResponse(path, media_type="application/gzip", filename="segmentation.nii.gz")
    except Exception as exc:
        logger.error(f"Error serving file: {exc}")
        raise HTTPException(status_code=500, detail="Error serving file") from exc


@router.get("/health")
async def health_check():
    return {"status": "healthy", "model_loaded": model is not None}
