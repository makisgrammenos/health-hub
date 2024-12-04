import torch
from PIL import Image
from .model.architecture import COVIDNext50
from .data.transforms import val_transforms
from .config import mapping, width, height

class COVIDNext50Model:
    def __init__(self, model_path: str):
        # Reverse mapping for predictions
        self.rev_mapping = {idx: name for name, idx in mapping.items()}
        
        # Load the model
        self.model = COVIDNext50(n_classes=len(self.rev_mapping))
        weights = torch.load(model_path, map_location=torch.device('cpu'))['state_dict']
        self.model.load_state_dict(weights)
        self.model.eval()

        # Define transforms
        self.transforms = val_transforms(width=width, height=height)

    def predict(self, img: Image.Image):
        # Apply transforms
        img_tensor = self.transforms(img).unsqueeze(0)

        # Perform prediction
        with torch.no_grad():
            logits = self.model(img_tensor)
            cat_id = int(torch.argmax(logits))
        return self.rev_mapping[cat_id]
