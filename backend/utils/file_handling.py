import zipfile
from pathlib import Path
from PIL import Image
from io import BytesIO
from typing import Dict
from fastapi import UploadFile
from pathlib import Path

def extract_zip(file):
    temp_dir = Path("temp_batch")
    temp_dir.mkdir(exist_ok=True)
    with zipfile.ZipFile(BytesIO(file.file.read()), "r") as zip_ref:
        zip_ref.extractall(temp_dir)
    return temp_dir

def load_image(file):
    return Image.open(file.file).convert("RGB")
import os

def save_uploaded_files(files: Dict[str, UploadFile], temp_dir: str) -> Dict[str, str]:
    """
    Save uploaded files to a temporary directory and return their paths.
    """
    file_paths = {}
    for modality, file in files.items():
        file_path = Path(temp_dir) / file.filename
        with open(file_path, "wb") as f:
            f.write(file.file.read())
        file_paths[modality] = str(file_path)
    return file_paths
