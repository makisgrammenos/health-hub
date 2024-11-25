import os
from typing import Dict
from fastapi import UploadFile
from pathlib import Path

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
6