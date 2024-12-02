from fastapi import APIRouter, WebSocket, WebSocketDisconnect
import cv2
import numpy as np
from PIL import Image
from io import BytesIO
import base64
import json

router = APIRouter()

# Utility functions
def image_to_bytes(image: np.ndarray):
    _, buffer = cv2.imencode('.png', image)
    return buffer.tobytes()

def process_image(image: np.ndarray, params: dict, methods: dict):
    # Validate and adjust parameters
    def ensure_odd(value):
        return int(value) if int(value) % 2 == 1 else int(value) + 1

    bilateral_d = ensure_odd(params.get("bilateralD", 9))
    bilateral_sigma_color = params.get("bilateralSigmaColor", 75)
    bilateral_sigma_space = params.get("bilateralSigmaSpace", 75)
    gaussian_kernel_size = ensure_odd(params.get("gaussianKernelSize", 5))
    clip_limit = params.get("clipLimit", 2.0)
    tile_grid_size = params.get("tileGridSize", 8)
    unsharp_strength = params.get("unsharpStrength", 1.5)
    normalize_min = params.get("normalizeMin", 0)
    normalize_max = params.get("normalizeMax", 255)

    # Processing Steps

    # Step 1: Denoising
    if methods.get('denoising', True):
        denoise_method = params.get("denoiseMethod", "bilateral")
        if denoise_method == "bilateral":
            image = cv2.bilateralFilter(
                image,
                d=bilateral_d,
                sigmaColor=bilateral_sigma_color,
                sigmaSpace=bilateral_sigma_space
            )
        elif denoise_method == "gaussian":
            image = cv2.GaussianBlur(image, (gaussian_kernel_size, gaussian_kernel_size), 0)
        elif denoise_method == "median":
            image = cv2.medianBlur(image, gaussian_kernel_size)

    # Step 2: Normalize Intensity
    if methods.get('normalization', True):
        image = cv2.normalize(image, None, alpha=normalize_min, beta=normalize_max, norm_type=cv2.NORM_MINMAX)

    # Step 3: Histogram Equalization
    if methods.get('histogramEqualization', True):
        image = cv2.equalizeHist(image)

    # Step 4: Contrast Enhancement (CLAHE)
    if methods.get('clahe', True):
        clahe = cv2.createCLAHE(clipLimit=clip_limit, tileGridSize=(tile_grid_size, tile_grid_size))
        image = clahe.apply(image)

    # Step 5: Edge Enhancement (Unsharp Masking)
    if methods.get('unsharpMasking', True):
        blurred = cv2.GaussianBlur(image, (gaussian_kernel_size, gaussian_kernel_size), 0)
        image = cv2.addWeighted(image, 1 + unsharp_strength, blurred, -unsharp_strength, 0)

    # Step 6: High-Pass Filtering (Optional)
    # if methods.get('highPassFiltering', False):
    #     highpass_kernel = np.array([[-1, -1, -1], [-1, 8, -1], [-1, -1, -1]], dtype=np.float32)
    #     image = cv2.filter2D(image, -1, highpass_kernel)

    return image

# WebSocket route
@router.websocket("/ws")
async def websocket_endpoint(websocket: WebSocket):
    await websocket.accept()

    try:
        while True:
            # Receive image and parameters from the client
            message = await websocket.receive_text()
            message_data = json.loads(message)

            # Parse the image data and parameters
            image_data = message_data.get("image")
            params = message_data.get("params", {})
            methods = message_data.get("methods", {})

            # Correctly decode the base64 image data
            image_bytes = base64.b64decode(image_data)
            pil_image = Image.open(BytesIO(image_bytes)).convert("L")
            image = np.array(pil_image)

            # Process the image with parameters and methods
            processed_image = process_image(image, params, methods)
            image_bytes = image_to_bytes(processed_image)
            await websocket.send_bytes(image_bytes)

            # Notify client that processing is complete
            await websocket.send_text(json.dumps({"message": "Processing complete"}))
    except WebSocketDisconnect:
        print("WebSocket disconnected")
    except Exception as e:
        await websocket.send_text(json.dumps({"error": str(e)}))
