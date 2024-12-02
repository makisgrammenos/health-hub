from fastapi import FastAPI, WebSocket, WebSocketDisconnect, APIRouter
from fastapi.middleware.cors import CORSMiddleware
import cv2
import numpy as np
from PIL import Image
from io import BytesIO
import base64
import json

router = APIRouter()

# Utility functions
def image_to_bytes(image: np.ndarray) -> bytes:
    """Convert an RGB image to PNG bytes."""
    _, buffer = cv2.imencode(".png", cv2.cvtColor(image, cv2.COLOR_RGB2BGR))
    return buffer.tobytes()

def apply_to_channels(image: np.ndarray, func):
    """Apply a function channel-wise to an RGB image."""
    channels = cv2.split(image)
    processed_channels = [func(channel) for channel in channels]
    return cv2.merge(processed_channels)

def process_image(image: np.ndarray, params: dict, methods: dict) -> np.ndarray:
    """Process the image based on user-selected methods."""
    def ensure_odd(value):
        return int(value) if int(value) % 2 == 1 else int(value) + 1

    # Extract parameters
    bilateral_d = ensure_odd(params.get("bilateralD", 9))
    bilateral_sigma_color = params.get("bilateralSigmaColor", 75)
    bilateral_sigma_space = params.get("bilateralSigmaSpace", 75)
    gaussian_kernel_size = ensure_odd(params.get("gaussianKernelSize", 5))
    clip_limit = params.get("clipLimit", 2.0)
    tile_grid_size = params.get("tileGridSize", 8)
    unsharp_strength = params.get("unsharpStrength", 1.5)
    normalize_min = params.get("normalizeMin", 0)
    normalize_max = params.get("normalizeMax", 255)
    window_width = params.get("windowWidth", 255)
    window_center = params.get("windowCenter", 127)
    threshold_min = params.get("thresholdMin", 50)
    threshold_max = params.get("thresholdMax", 200)
    pseudocolor_map = getattr(cv2, params.get("pseudocolorMap", "COLORMAP_JET"), cv2.COLORMAP_JET)
    kernel_size = ensure_odd(params.get("kernelSize", 3))
    morph_operation = params.get("morphOperation", "dilation")

    # Step 1: Denoising
    if methods.get("denoising", False):
        denoise_method = params.get("denoiseMethod", "bilateral")
        if denoise_method == "bilateral":
            image = cv2.bilateralFilter(image, bilateral_d, bilateral_sigma_color, bilateral_sigma_space)
        elif denoise_method == "gaussian":
            image = cv2.GaussianBlur(image, (gaussian_kernel_size, gaussian_kernel_size), 0)
        elif denoise_method == "median":
            image = cv2.medianBlur(image, gaussian_kernel_size)

    # Step 2: Normalize Intensity
    if methods.get("normalization", False):
        image = cv2.normalize(image, None, alpha=normalize_min, beta=normalize_max, norm_type=cv2.NORM_MINMAX)

    # Step 3: Histogram Equalization
    if methods.get("histogramEqualization", False):
        image = apply_to_channels(image, cv2.equalizeHist)

    # Step 4: Contrast Enhancement (CLAHE)
    if methods.get("clahe", False):
        clahe = cv2.createCLAHE(clipLimit=clip_limit, tileGridSize=(tile_grid_size, tile_grid_size))
        image = apply_to_channels(image, clahe.apply)

    # Step 5: Edge Enhancement (Unsharp Masking)
    if methods.get("unsharpMasking", False):
        blurred = cv2.GaussianBlur(image, (gaussian_kernel_size, gaussian_kernel_size), 0)
        image = cv2.addWeighted(image, 1 + unsharp_strength, blurred, -unsharp_strength, 0)

    # Step 6: Windowing and Leveling
    if methods.get("windowing", False):
        min_val = window_center - (window_width // 2)
        max_val = window_center + (window_width // 2)
        image = np.clip(image, min_val, max_val)
        image = cv2.normalize(image, None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)

    # Step 7: Threshold Segmentation
    if methods.get("thresholdSegmentation", False):
        _, image = cv2.threshold(image, threshold_min, threshold_max, cv2.THRESH_BINARY)

    # Step 8: Pseudocoloring
    if methods.get("pseudocolor", False):
        image = cv2.applyColorMap(image, pseudocolor_map)

    # Step 9: Morphological Operations
    if methods.get("morphologicalOperations", False):
        kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (kernel_size, kernel_size))
        if morph_operation == "dilation":
            image = cv2.dilate(image, kernel)
        elif morph_operation == "erosion":
            image = cv2.erode(image, kernel)
        elif morph_operation == "opening":
            image = cv2.morphologyEx(image, cv2.MORPH_OPEN, kernel)
        elif morph_operation == "closing":
            image = cv2.morphologyEx(image, cv2.MORPH_CLOSE, kernel)

    return image

@router.websocket("/ws")
async def websocket_endpoint(websocket: WebSocket):
    await websocket.accept()
    try:
        while True:
            # Receive data
            message = await websocket.receive_text()
            data = json.loads(message)
            image_data = data.get("image")
            params = data.get("params", {})
            methods = data.get("methods", {})

            # Decode image
            image_bytes = base64.b64decode(image_data)
            pil_image = Image.open(BytesIO(image_bytes)).convert("RGB")
            image = np.array(pil_image)

            # Process image
            processed_image = process_image(image, params, methods)

            # Send processed image back
            processed_bytes = image_to_bytes(processed_image)
            base64_image = base64.b64encode(processed_bytes).decode()
            await websocket.send_text(base64_image)
    except WebSocketDisconnect:
        print("WebSocket disconnected")
    except Exception as e:
        await websocket.send_text(json.dumps({"error": str(e)}))
