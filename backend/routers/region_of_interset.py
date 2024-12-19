from fastapi import APIRouter, WebSocket, WebSocketDisconnect, HTTPException
from fastapi.responses import FileResponse
from PIL import Image
import io
import base64
import os
import uuid
from pathlib import Path

CROPPED_IMAGE_DIR = "cropped_images"
os.makedirs(CROPPED_IMAGE_DIR, exist_ok=True)

router = APIRouter()

@router.websocket("/ws")
async def websocket_endpoint(websocket: WebSocket):
    await websocket.accept()
    try:
        while True:
            data = await websocket.receive_json()
            action = data.get("action")

            if action == "upload_image":
                image_base64 = data.get("image")
                try:
                    image_data = base64.b64decode(image_base64)
                    image = Image.open(io.BytesIO(image_data)).convert("RGB")
                    # Generate a unique filename to prevent overwriting
                    unique_filename = f"uploaded_{uuid.uuid4().hex}.jpg"
                    image.save(unique_filename)
                    await websocket.send_json({"message": "Image uploaded successfully", "filename": unique_filename})
                except Exception as e:
                    await websocket.send_json({"error": f"Failed to upload image: {str(e)}"})

            elif action == "crop":
                roi = data.get("roi")
                filename = data.get("filename")  # Receive the filename of the uploaded image
                if not filename:
                    await websocket.send_json({"error": "No filename provided for cropping"})
                    continue

                x, y, width, height = int(roi["x"]), int(roi["y"]), int(roi["width"]), int(roi["height"])

                uploaded_image_path = Path(filename)
                if not uploaded_image_path.exists():
                    await websocket.send_json({"error": "No uploaded image found with the provided filename"})
                    continue

                try:
                    # Open the uploaded image
                    image = Image.open(uploaded_image_path)
                    cropped_image = image.crop((x, y, x + width, y + height))

                    # Generate a unique filename for the cropped image
                    cropped_filename = f"cropped_{uuid.uuid4().hex}.jpg"
                    cropped_path = Path(CROPPED_IMAGE_DIR) / cropped_filename
                    cropped_image.save(cropped_path)

                    # Encode cropped image as base64
                    img_byte_arr = io.BytesIO()
                    cropped_image.save(img_byte_arr, format='JPEG')
                    img_base64 = base64.b64encode(img_byte_arr.getvalue()).decode('utf-8')

                    await websocket.send_json({
                        "message": "Image cropped successfully",
                        "image_base64": img_base64,
                        "cropped_filename": cropped_filename
                    })
                except Exception as e:
                    await websocket.send_json({"error": f"Failed to crop image: {str(e)}"})

    except WebSocketDisconnect:
        print("Client disconnected")
    except Exception as e:
        # Catch any other exceptions to prevent the server from crashing
        print(f"Unexpected error: {e}")
        await websocket.close()

@router.get("/download/{filename}")
async def download_cropped_image(filename: str):
    file_path = Path(CROPPED_IMAGE_DIR) / filename
    if not file_path.exists():
        raise HTTPException(status_code=404, detail="File not found")

    return FileResponse(path=file_path, media_type='image/jpeg', filename=filename)

# New Endpoint to Delete Cropped Images
@router.delete("/delete/{filename}")
async def delete_cropped_image(filename: str):
    file_path = Path(CROPPED_IMAGE_DIR) / filename
    if not file_path.exists():
        raise HTTPException(status_code=404, detail="File not found")
    
    try:
        os.remove(file_path)
        return {"message": f"File '{filename}' deleted successfully."}
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to delete file: {str(e)}")
