import nibabel as nib
import numpy as np


def save_nifti(prediction: np.ndarray, output_path: str, affine=None):
    """
    Save the segmentation prediction as a NIfTI file.

    Args:
        prediction: The segmentation prediction array
        output_path: Path where to save the NIfTI file
        affine: Affine transformation matrix for the NIfTI image (default: identity matrix)
    """
    if affine is None:
        affine = np.eye(4)

    # Ensure prediction is in correct format (labels should be [0, 1, 2] for background, tumor core, edema)
    # For binary masks, we need to ensure values are 0 or 1
    prediction = prediction.astype(np.uint8)

    nifti_image = nib.Nifti1Image(prediction, affine=affine)
    nib.save(nifti_image, output_path)