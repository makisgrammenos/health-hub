from fastapi import APIRouter, UploadFile, File, HTTPException, BackgroundTasks
from fastapi.responses import FileResponse
from typing import List, Optional
import tempfile
import torch
import os
import numpy as np
import base64
import io
import logging
from PIL import Image
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

# MONAI imports
from monai.transforms import (
    Compose,
    LoadImaged,
    EnsureChannelFirstd,
    ScaleIntensityd,
    EnsureTyped,
)
from monai.data import PydicomReader, ITKReader
import nibabel as nib
import pydicom

# Import your model class
from models.varscular_seg.model import MonaiUNet2DSegmentationModel
from models.varscular_seg.config import MODEL_CHECKPOINT_PATH

# --------------------------------------------------------------------------- #
# Logging & Global Initialization
# --------------------------------------------------------------------------- #
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

router = APIRouter()
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# Initialize model
try:
    model = MonaiUNet2DSegmentationModel(MODEL_CHECKPOINT_PATH, device)
    logger.info(f"MONAI UNet model loaded successfully on {device}")
except Exception as e:
    logger.error(f"Failed to load MONAI UNet model: {e},")
    model = None

# Colormap for segmentation visualization (4 classes)
seg_colors = [
    (0, 0, 0, 0),        # 0: background
    (1, 0, 0, 0.7),      # 1: class 1 (red)
    (0, 1, 0, 0.7),      # 2: class 2 (green)
    (0, 0, 1, 0.7),      # 3: class 3 (blue)
]
seg_cmap = ListedColormap(seg_colors, name="seg_cmap")
seg_norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], ncolors=len(seg_colors))

# --------------------------------------------------------------------------- #
# Helper Functions
# --------------------------------------------------------------------------- #
def get_dicom_reader():
    """Get appropriate DICOM reader"""
    try:
        return PydicomReader()
    except ImportError:
        try:
            return ITKReader()
        except ImportError:
            logger.warning("No specialized DICOM reader available, using default")
            return None


def create_overlay_image(image_slice: np.ndarray, seg_slice: np.ndarray) -> str:
    """
    Create base64-encoded PNG with segmentation overlay
    
    Args:
        image_slice: 2D grayscale image
        seg_slice: 2D segmentation mask
    
    Returns:
        Base64-encoded PNG string
    """
    # Normalize image
    img_norm = (
        (image_slice - image_slice.min()) / (image_slice.ptp() + 1e-8)
        if image_slice.max() > image_slice.min()
        else np.zeros_like(image_slice)
    )
    
    # Create figure
    plt.figure(figsize=(10, 10), dpi=100)
    plt.imshow(img_norm, cmap="gray")
    plt.imshow(seg_slice, cmap=seg_cmap, norm=seg_norm, alpha=0.5)
    plt.axis("off")
    plt.tight_layout(pad=0)
    
    # Convert to base64
    buf = io.BytesIO()
    plt.savefig(buf, format="png", bbox_inches="tight", pad_inches=0)
    plt.close()
    buf.seek(0)
    
    return "data:image/png;base64," + base64.b64encode(buf.read()).decode("utf-8")


def save_nifti(prediction: np.ndarray, output_path: str, affine=None):
    """Save segmentation as NIfTI file"""
    if affine is None:
        affine = np.eye(4)
    
    prediction = prediction.astype(np.uint8)
    nifti_image = nib.Nifti1Image(prediction, affine=affine)
    nib.save(nifti_image, output_path)


def process_dicom_file(file_path: str) -> tuple:
    """
    Process a single DICOM file
    
    Returns:
        Tuple of (image_array, metadata)
    """
    # Try pydicom first
    try:
        ds = pydicom.dcmread(file_path)
        image = ds.pixel_array.astype(np.float32)
        
        # Apply rescale if available
        if hasattr(ds, 'RescaleSlope') and hasattr(ds, 'RescaleIntercept'):
            image = image * ds.RescaleSlope + ds.RescaleIntercept
        
        metadata = {
            'PatientID': getattr(ds, 'PatientID', 'Unknown'),
            'StudyDate': getattr(ds, 'StudyDate', 'Unknown'),
            'Modality': getattr(ds, 'Modality', 'Unknown'),
            'SliceLocation': getattr(ds, 'SliceLocation', 0),
            'InstanceNumber': getattr(ds, 'InstanceNumber', 0),
        }
        
        return image, metadata
    except Exception as e:
        logger.error(f"Error reading DICOM with pydicom: {e}")
        raise


# --------------------------------------------------------------------------- #
# API Routes
# --------------------------------------------------------------------------- #
# @router.post("/segment")
# async def segment_dicom(
#     background_tasks: BackgroundTasks,
#     files: List[UploadFile] = File(..., description="DICOM files to segment"),
#     return_all_slices: bool = True
# ):
#     """
#     Segment DICOM images using MONAI UNet
    
#     Args:
#         files: List of DICOM files
#         return_all_slices: If True, return all slices; if False, return representative slices
#     """
#     if model is None:
#         raise HTTPException(
#             status_code=503,
#             detail="Model is not available. Please try again later."
#         )
    
#     temp_files = []
#     output_paths = []
    
#     try:
#         # ------------------------------------------------------------------ #
#         # 1. Save uploaded files and prepare for processing
#         # ------------------------------------------------------------------ #
#         dicom_data = []
        
#         for file in files:
#             # Validate file extension
#             if not file.filename.lower().endswith(('.dcm', '.dicom')):
#                 # Check if it might be a DICOM without extension
#                 logger.warning(f"File {file.filename} doesn't have standard DICOM extension")
            
#             # Save to temporary file
#             with tempfile.NamedTemporaryFile(delete=False, suffix='.dcm') as tmp:
#                 content = await file.read()
#                 tmp.write(content)
#                 temp_files.append(tmp.name)
                
#                 # Process DICOM
#                 try:
#                     image, metadata = process_dicom_file(tmp.name)
#                     dicom_data.append({
#                         'image': image,
#                         'metadata': metadata,
#                         'filename': file.filename
#                     })
#                 except Exception as e:
#                     logger.error(f"Failed to process {file.filename}: {e}")
#                     continue
        
#         if not dicom_data:
#             raise HTTPException(
#                 status_code=400,
#                 detail="No valid DICOM files could be processed"
#             )
        
#         logger.info(f"Processing {len(dicom_data)} DICOM file(s)")
        
#         # ------------------------------------------------------------------ #
#         # 2. Sort by slice location if multiple files
#         # ------------------------------------------------------------------ #
#         if len(dicom_data) > 1:
#             dicom_data.sort(key=lambda x: x['metadata'].get('SliceLocation', 0))
        
#         # ------------------------------------------------------------------ #
#         # 3. Process each slice
#         # ------------------------------------------------------------------ #
#         segmentation_results = []
#         all_segmentations = []
        
#         for idx, data in enumerate(dicom_data):
#             image = data['image']
            
#             # Ensure image is 2D
#             if image.ndim == 3:
#                 image = image[..., 0]  # Take first channel if multi-channel
            
#             # Normalize image
#             image = (image - image.min()) / (image.ptp() + 1e-8) if image.ptp() > 0 else image
            
#             # Convert to tensor and add dimensions (C, H, W)
#             image_tensor = torch.from_numpy(image).float().unsqueeze(0)
            
#             # Run prediction
#             seg_pred = model.predict_single_slice(image_tensor)
            
#             # Store segmentation
#             all_segmentations.append(seg_pred)
            
#             # Create visualization
#             if return_all_slices or idx % max(1, len(dicom_data) // 3) == 0:
#                 overlay = create_overlay_image(image, seg_pred)
#                 segmentation_results.append({
#                     'slice_index': idx,
#                     'filename': data['filename'],
#                     'image': overlay
#                 })
        
#         # ------------------------------------------------------------------ #
#         # 4. Create 3D volume if multiple slices
#         # ------------------------------------------------------------------ #
#         if len(all_segmentations) > 1:
#             # Stack into 3D volume
#             volume_seg = np.stack(all_segmentations, axis=-1)
#         else:
#             volume_seg = all_segmentations[0][..., np.newaxis]
        
#         # Save as NIfTI
#         output_file = tempfile.NamedTemporaryFile(delete=False, suffix='.nii.gz')
#         output_path = output_file.name
#         output_file.close()
#         save_nifti(volume_seg, output_path)
#         output_paths.append(output_path)
        
#         # ------------------------------------------------------------------ #
#         # 5. Calculate statistics
#         # ------------------------------------------------------------------ #
#         unique_labels, counts = np.unique(volume_seg, return_counts=True)
#         label_stats = {}
        
#         for label, count in zip(unique_labels, counts):
#             if label == 0:
#                 label_stats['background'] = int(count)
#             else:
#                 label_stats[f'class_{int(label)}'] = int(count)
        
#         total_voxels = volume_seg.size
#         label_percentages = {
#             k: f"{(v/total_voxels)*100:.2f}%" for k, v in label_stats.items()
#         }
        
#         # ------------------------------------------------------------------ #
#         # 6. Schedule cleanup and return results
#         # ------------------------------------------------------------------ #
#         for path in temp_files:
#             background_tasks.add_task(os.unlink, path)
#         for path in output_paths:
#             background_tasks.add_task(os.unlink, path)
        
#         return {
#             'status': 'success',
#             'segmented_slices': segmentation_results,
#             'download_url': f'/imaging/monai-unet/download/{os.path.basename(output_path)}',
#             'statistics': {
#                 'total_slices': len(dicom_data),
#                 'volume_shape': list(volume_seg.shape),
#                 'label_counts': label_stats,
#                 'label_percentages': label_percentages
#             }
#         }
        
#     except Exception as e:
#         logger.error(f"Error during segmentation: {e}", exc_info=True)
        
#         # Cleanup
#         for path in temp_files:
#             if os.path.exists(path):
#                 try:
#                     os.unlink(path)
#                 except:
#                     pass
#         for path in output_paths:
#             if os.path.exists(path):
#                 try:
#                     os.unlink(path)
#                 except:
#                     pass
        
#         raise HTTPException(
#             status_code=500,
#             detail=f"Error during segmentation: {str(e)}"
#         )
import zipfile
import tempfile
import shutil

@router.post("/segment")
async def segment_dicom(
    background_tasks: BackgroundTasks,
    files: List[UploadFile] = File(..., description="DICOM files or a ZIP containing DICOMs"),
    return_all_slices: bool = True
):
    """
    Segment DICOM images using MONAI UNet.
    Accepts either:
      - Multiple .dcm/.dicom files directly, OR
      - A single .zip containing multiple DICOM slices
    """
    if model is None:
        raise HTTPException(status_code=503, detail="Model is not available. Please try again later.")

    temp_files = []
    output_paths = []
    extracted_dir = None

    try:
        dicom_files = []

        # ------------------------------------------------------------------ #
        # 1. Check if user uploaded a ZIP
        # ------------------------------------------------------------------ #
        if len(files) == 1 and files[0].filename.lower().endswith('.zip'):
            # Create a temp dir for extracted files
            extracted_dir = tempfile.mkdtemp()

            # Save uploaded ZIP
            zip_path = os.path.join(extracted_dir, "upload.zip")
            with open(zip_path, "wb") as zf:
                content = await files[0].read()
                zf.write(content)

            # Extract ZIP
            with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                zip_ref.extractall(extracted_dir)

            # Collect all .dcm/.dicom files from extracted dir
            for root, _, filenames in os.walk(extracted_dir):
                for fname in filenames:
                    if fname.lower().endswith((".dcm", ".dicom")):
                        dicom_files.append(os.path.join(root, fname))

            if not dicom_files:
                raise HTTPException(status_code=400, detail="ZIP does not contain any valid DICOM files.")

        else:
            # ------------------------------------------------------------------ #
            # 2. Handle normal DICOM upload
            # ------------------------------------------------------------------ #
            for file in files:
                if not file.filename.lower().endswith((".dcm", ".dicom")):
                    logger.warning(f"File {file.filename} does not have a standard DICOM extension")
                with tempfile.NamedTemporaryFile(delete=False, suffix='.dcm') as tmp:
                    content = await file.read()
                    tmp.write(content)
                    temp_files.append(tmp.name)
                    dicom_files.append(tmp.name)

        if not dicom_files:
            raise HTTPException(status_code=400, detail="No valid DICOM files to process.")

        logger.info(f"Processing {len(dicom_files)} DICOM file(s)")

        # ------------------------------------------------------------------ #
        # 3. Process DICOM files
        # ------------------------------------------------------------------ #
        dicom_data = []
        for path in dicom_files:
            try:
                image, metadata = process_dicom_file(path)
                dicom_data.append({
                    'image': image,
                    'metadata': metadata,
                    'filename': os.path.basename(path)
                })
            except Exception as e:
                logger.error(f"Failed to process {path}: {e}")
                continue

        if not dicom_data:
            raise HTTPException(status_code=400, detail="No valid DICOM data could be processed.")

        # Sort slices if multiple
        if len(dicom_data) > 1:
            dicom_data.sort(key=lambda x: x['metadata'].get('SliceLocation', 0))

        segmentation_results = []
        all_segmentations = []

        for idx, data in enumerate(dicom_data):
            image = data['image']
            if image.ndim == 3:
                image = image[..., 0]
            image = (image - image.min()) / (image.ptp() + 1e-8) if image.ptp() > 0 else image
            image_tensor = torch.from_numpy(image).float().unsqueeze(0)
            seg_pred = model.predict_single_slice(image_tensor)
            all_segmentations.append(seg_pred)

            if return_all_slices or idx % max(1, len(dicom_data)//3) == 0:
                overlay = create_overlay_image(image, seg_pred)
                segmentation_results.append({
                    'slice_index': idx,
                    'filename': data['filename'],
                    'image': overlay
                })

        # ------------------------------------------------------------------ #
        # 4. Save NIfTI result
        # ------------------------------------------------------------------ #
        volume_seg = np.stack(all_segmentations, axis=-1) if len(all_segmentations) > 1 else all_segmentations[0][..., np.newaxis]
        output_file = tempfile.NamedTemporaryFile(delete=False, suffix='.nii.gz')
        output_path = output_file.name
        output_file.close()
        save_nifti(volume_seg, output_path)
        output_paths.append(output_path)

        # ------------------------------------------------------------------ #
        # 5. Statistics
        # ------------------------------------------------------------------ #
        unique_labels, counts = np.unique(volume_seg, return_counts=True)
        label_stats = {('background' if l == 0 else f'class_{int(l)}'): int(c) for l, c in zip(unique_labels, counts)}
        total_voxels = volume_seg.size
        label_percentages = {k: f"{(v/total_voxels)*100:.2f}%" for k, v in label_stats.items()}

        # ------------------------------------------------------------------ #
        # 6. Cleanup
        # ------------------------------------------------------------------ #
        for path in temp_files:
            background_tasks.add_task(os.unlink, path)
        for path in output_paths:
            background_tasks.add_task(os.unlink, path)
        if extracted_dir:
            background_tasks.add_task(shutil.rmtree, extracted_dir)

        return {
            'status': 'success',
            'segmented_slices': segmentation_results,
            'download_url': f'/imaging/monai-unet/download/{os.path.basename(output_path)}',
            'statistics': {
                'total_slices': len(dicom_data),
                'volume_shape': list(volume_seg.shape),
                'label_counts': label_stats,
                'label_percentages': label_percentages
            }
        }

    except Exception as e:
        logger.error(f"Error during segmentation: {e}", exc_info=True)
        if extracted_dir and os.path.exists(extracted_dir):
            shutil.rmtree(extracted_dir)
        raise HTTPException(status_code=500, detail=f"Error during segmentation: {str(e)}")

    
@router.post("/segment-nifti")
async def segment_nifti(
    background_tasks: BackgroundTasks,
    file: UploadFile = File(..., description="NIfTI file to segment")
):
    """
    Segment a NIfTI volume using MONAI UNet
    """
    if model is None:
        raise HTTPException(
            status_code=503,
            detail="Model is not available. Please try again later."
        )
    
    temp_file = None
    output_path = None
    
    try:
        # Validate file extension
        if not file.filename.lower().endswith(('.nii', '.nii.gz')):
            raise HTTPException(
                status_code=400,
                detail="Invalid file format. Only .nii or .nii.gz files are supported."
            )
        
        # Save uploaded file
        suffix = '.nii.gz' if file.filename.endswith('.gz') else '.nii'
        with tempfile.NamedTemporaryFile(delete=False, suffix=suffix) as tmp:
            content = await file.read()
            tmp.write(content)
            temp_file = tmp.name
        
        # Load NIfTI
        nifti_img = nib.load(temp_file)
        volume = nifti_img.get_fdata()
        affine = nifti_img.affine
        
        # Process volume
        if volume.ndim == 2:
            volume = volume[..., np.newaxis]  # Add depth dimension
        
        # Normalize
        volume = (volume - volume.min()) / (volume.ptp() + 1e-8) if volume.ptp() > 0 else volume
        
        # Convert to tensor (C, D, H, W) for 3D processing
        volume_tensor = torch.from_numpy(volume).float()
        volume_tensor = volume_tensor.permute(2, 0, 1).unsqueeze(0).unsqueeze(0)  # (B, C, D, H, W)
        
        # Run prediction
        seg_pred = model.predict(volume_tensor)
        
        # Reorder dimensions back to (H, W, D)
        if seg_pred.ndim == 3:
            seg_pred = seg_pred.transpose(1, 2, 0)
        
        # Save segmentation
        output_file = tempfile.NamedTemporaryFile(delete=False, suffix='.nii.gz')
        output_path = output_file.name
        output_file.close()
        save_nifti(seg_pred, output_path, affine=affine)
        
        # Generate preview slices
        preview_slices = []
        slice_indices = [0, seg_pred.shape[2]//2, seg_pred.shape[2]-1]
        
        for idx in slice_indices:
            if 0 <= idx < seg_pred.shape[2]:
                overlay = create_overlay_image(
                    volume[:, :, idx],
                    seg_pred[:, :, idx]
                )
                preview_slices.append({
                    'slice_index': idx,
                    'image': overlay
                })
        
        # Calculate statistics
        unique_labels, counts = np.unique(seg_pred, return_counts=True)
        label_stats = {f'class_{int(label)}': int(count) for label, count in zip(unique_labels, counts)}
        
        # Cleanup
        background_tasks.add_task(os.unlink, temp_file)
        background_tasks.add_task(os.unlink, output_path)
        
        return {
            'status': 'success',
            'preview_slices': preview_slices,
            'download_url': f'/imaging/monai-unet/download/{os.path.basename(output_path)}',
            'statistics': {
                'volume_shape': list(seg_pred.shape),
                'label_counts': label_stats
            }
        }
        
    except Exception as e:
        logger.error(f"Error during NIfTI segmentation: {e}", exc_info=True)
        
        # Cleanup
        if temp_file and os.path.exists(temp_file):
            try:
                os.unlink(temp_file)
            except:
                pass
        if output_path and os.path.exists(output_path):
            try:
                os.unlink(output_path)
            except:
                pass
        
        raise HTTPException(
            status_code=500,
            detail=f"Error during segmentation: {str(e)}"
        )


@router.get("/download/{filename}")
async def download_segmentation(filename: str):
    """Download segmentation result"""
    path = os.path.join(tempfile.gettempdir(), filename)
    
    if not os.path.exists(path):
        raise HTTPException(
            status_code=404,
            detail="File not found or has expired"
        )
    
    try:
        return FileResponse(
            path,
            media_type="application/gzip",
            filename=f"segmentation_{filename}"
        )
    except Exception as e:
        logger.error(f"Error serving file: {e}")
        raise HTTPException(
            status_code=500,
            detail="Error serving file"
        )


@router.get("/health")
async def health_check():
    """Health check endpoint"""
    return {
        "status": "healthy",
        "model_loaded": model is not None,
        "device": str(device)
    }