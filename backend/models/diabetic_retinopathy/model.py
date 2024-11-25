import torch
from .preprocessing import preprocess_image
from .config import MODEL_PATH, CLASS_NAMES
# from .response import build_response
import numpy as np
import  torch.nn as nn
import  torchvision
from utils.file_handling import load_image

# Load the model
class InceptionModel(nn.Module):
    def __init__(self, transfer=False, *args, **kwargs):

        super(InceptionModel, self).__init__()
        self.model = torchvision.models.inception_v3(weights='DEFAULT')
        self.model.aux_logits = False
        if transfer:
            for param in self.model.parameters():
                param.requires_grad = False

        in_features = self.model.fc.in_features

        self.model.fc = nn.Linear(in_features, 2)

    def forward(self, x):

        x = self.model(x)
        return x


def load_model():
    checkpoint = torch.load(MODEL_PATH, map_location=torch.device("cpu"))
    model = InceptionModel()
    model.load_state_dict(checkpoint['model_state_dict'])
    model.eval()
    return model

model = load_model()

async def predict_single(file):
    image = load_image(file)
    preprocessed_image = preprocess_image(image)
    with torch.no_grad():
        outputs = model(preprocessed_image)
        _, predicted_class = torch.max(outputs, 1)
        confidence = torch.softmax(outputs, dim=1)[0, predicted_class].item()
        predicted_label = CLASS_NAMES[predicted_class.item()]
    prediction = {"result": predicted_label, "confidence": confidence}
    return prediction


async def predict_batch(temp_dir):
    predictions = []
    for image_path in temp_dir.iterdir():
        if image_path.suffix in [".jpeg", ".jpg", ".png"]:
            image = load_image(image_path)
            preprocessed_image = preprocess_image(image)
            with torch.no_grad():
                outputs = model(preprocessed_image)
                _, predicted_class = torch.max(outputs, 1)
                confidence = torch.softmax(outputs, dim=1)[0, predicted_class].item()
                predicted_label = CLASS_NAMES[predicted_class.item()]
            predictions.append({"image": image_path.name, "label": predicted_label, "confidence": confidence})
    return predictions
