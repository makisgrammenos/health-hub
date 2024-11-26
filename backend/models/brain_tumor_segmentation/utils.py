import nibabel as nib
import numpy as np


def save_nifti(prediction: np.ndarray, output_path: str):
    nifti_image = nib.Nifti1Image(prediction.astype(np.uint8), affine=np.eye(4))
    nib.save(nifti_image, output_path)
