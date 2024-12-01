from fastapi import APIRouter, UploadFile, File, Form, HTTPException
from fastapi.responses import JSONResponse
import cv2
import numpy as np
from PIL import Image
from io import BytesIO
import base64
import json

router = APIRouter()

# Utility functions
def image_to_base64(image: np.ndarray):
    # Encode image to PNG format
    is_success, buffer = cv2.imencode('.png', image)
    if not is_success:
        raise ValueError("Could not encode image to bytes")
    base64_str = base64.b64encode(buffer).decode('utf-8')
    return base64_str

def process_image_in_steps(image: np.ndarray, params: dict):
    steps = []

    # Ensure image is uint8
    if image.dtype != np.uint8:
        image = image.astype(np.uint8)

    # Step 1: Denoising
    denoise_method = params.get("denoise_method", "bilateral")
    if denoise_method == "bilateral":
        image = cv2.bilateralFilter(
            image,
            d=int(params.get("bilateral_d", 9)),
            sigmaColor=float(params.get("bilateral_sigma_color", 75)),
            sigmaSpace=float(params.get("bilateral_sigma_space", 75))
        )
    elif denoise_method == "gaussian":
        kernel_size = int(params.get("gaussian_kernel_size", 5))
        kernel_size = kernel_size if kernel_size % 2 == 1 else kernel_size + 1
        image = cv2.GaussianBlur(image, (kernel_size, kernel_size), 0)
    elif denoise_method == "median":
        kernel_size = int(params.get("gaussian_kernel_size", 5))
        kernel_size = kernel_size if kernel_size % 2 == 1 else kernel_size + 1
        image = cv2.medianBlur(image, kernel_size)
    steps.append(image.copy())

    # Step 2: Normalize Intensity
    normalize_min = int(params.get("normalize_min", 0))
    normalize_max = int(params.get("normalize_max", 255))
    image = cv2.normalize(image, None, alpha=normalize_min, beta=normalize_max, norm_type=cv2.NORM_MINMAX)
    steps.append(image.copy())

    # Step 3: Histogram Equalization
    image = cv2.equalizeHist(image)
    steps.append(image.copy())

    # Step 4: Contrast Enhancement (CLAHE)
    clip_limit = float(params.get("clip_limit", 2.0))
    tile_grid_size = int(params.get("tile_grid_size", 8))
    clahe = cv2.createCLAHE(clipLimit=clip_limit, tileGridSize=(tile_grid_size, tile_grid_size))
    image = clahe.apply(image)
    steps.append(image.copy())

    # Step 5: Edge Enhancement (Unsharp Masking)
    gaussian_kernel_size = int(params.get("gaussian_kernel_size", 5))
    gaussian_kernel_size = gaussian_kernel_size if gaussian_kernel_size % 2 == 1 else gaussian_kernel_size + 1
    unsharp_strength = float(params.get("unsharp_strength", 1.5))
    blurred = cv2.GaussianBlur(image, (gaussian_kernel_size, gaussian_kernel_size), 0)
    image = cv2.addWeighted(image, 1 + unsharp_strength, blurred, -unsharp_strength, 0)
    steps.append(image.copy())

    # Step 6: High-Pass Filtering
    highpass_kernel = params.get("highpass_kernel", [[-1, -1, -1], [-1, 8, -1], [-1, -1, -1]])
    kernel = np.array(highpass_kernel, dtype=np.float32)
    image = cv2.filter2D(image, -1, kernel)
    steps.append(image.copy())

    return steps

# HTTP POST route
@router.post("/process-image")
async def process_image(
    file: UploadFile = File(...),
    params_json: str = Form(...)
):
    try:
        # Read the image file
        contents = await file.read()
        pil_image = Image.open(BytesIO(contents)).convert("L")
        image = np.array(pil_image)

        # Parse parameters
        params = json.loads(params_json)

        # For debugging: print received parameters
        print("Received parameters:", params)

        # Process the image step-by-step with parameters
        steps = process_image_in_steps(image, params)

        # Convert processed images to base64 strings
        processed_images = [image_to_base64(step) for step in steps]

        return JSONResponse(content={"images": processed_images})

    except Exception as e:
        error_message = f"Error processing image: {str(e)}"
        print(error_message)
        raise HTTPException(status_code=500, detail=error_message)
