# models/monai_unet_segmentation/config.py

import os

# Model configuration
MODEL_CHECKPOINT_PATH = os.environ.get(
    "MONAI_UNET_MODEL_PATH", 
    "models/varscular_seg/model.pt"
)

# Ensure the model path exists
if not os.path.exists(MODEL_CHECKPOINT_PATH):
    print(f"⚠️  Warning: Model checkpoint not found at {MODEL_CHECKPOINT_PATH}")
    print("Please ensure the model file is placed at the correct location.")

# Inference configuration
INFERENCE_CONFIG = {
    "roi_size": (256, 256),
    "sw_batch_size": 1,
    "overlap": 0.25,
    "num_classes": 4,
}

# Supported file formats
SUPPORTED_FORMATS = {
    "dicom": [".dcm", ".dicom"],
    "nifti": [".nii", ".nii.gz"]
}