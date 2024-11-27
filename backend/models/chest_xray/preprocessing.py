import torch
import numpy as np
import torchvision.transforms as transforms
import torchxrayvision as xrv

def preprocess_image(img, device):
    img = xrv.datasets.normalize(img, 255)  # Normalize to [-1024, 1024]
    if img.ndim == 3:
        img = img.mean(axis=2)  # Convert to grayscale
    img = img[None, ...]  # Add channel dimension
    transform = transforms.Compose([
        xrv.datasets.XRayCenterCrop(),
        xrv.datasets.XRayResizer(224)
    ])
    img = transform(img)
    img_tensor = torch.from_numpy(img).float().unsqueeze(0).to(device)  # Add batch dimension
    return img_tensor
