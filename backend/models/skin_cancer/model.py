import torch
from torchvision import models, transforms

def load_model(model_path: str, device: torch.device):
    """
    Load and initialize the model.

    Args:
        model_path (str): Path to the model's state_dict file.
        device (torch.device): Device to load the model onto.

    Returns:
        torch.nn.Module: Loaded and initialized model.
    """
    # Initialize the model
    model = models.efficientnet_v2_l(weights='DEFAULT')
    model.classifier[1] = torch.nn.Linear(model.classifier[1].in_features, 1)

    # Load the model state_dict
    state_dict = torch.load(model_path, map_location=device)
    model.load_state_dict(state_dict)
    model.eval()
    model.to(device)

    return model

def get_image_transform(img_size: int = 112):
    """
    Define the image transformation pipeline.

    Args:
        img_size (int): Target image size.

    Returns:
        torchvision.transforms.Compose: Transformation pipeline.
    """
    return transforms.Compose([
        transforms.Resize(size=(img_size, img_size), antialias=True),
        transforms.CenterCrop(size=(img_size, img_size)),
        transforms.ToTensor(),
        transforms.Normalize(mean=[0.485, 0.456, 0.406], std=[0.229, 0.224, 0.225])
    ])
