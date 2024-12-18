import torch
import torchvision.transforms as transforms
from PIL import Image

class BreastCancerModel:
    def __init__(self, model_path: str):
        # Initialize the model
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.model = torch.hub.load('pytorch/vision', 'resnet18', pretrained=False)
        num_features = self.model.fc.in_features
        self.model.fc = torch.nn.Sequential(
            torch.nn.Linear(num_features, 512),
            torch.nn.ReLU(),
            torch.nn.BatchNorm1d(512),
            torch.nn.Dropout(0.5),
            torch.nn.Linear(512, 256),
            torch.nn.ReLU(),
            torch.nn.BatchNorm1d(256),
            torch.nn.Dropout(0.5),
            torch.nn.Linear(256, 2)
        )
        self.model.load_state_dict(torch.load(model_path, map_location=self.device))
        self.model = self.model.to(self.device)
        self.model.eval()

        # Define preprocessing transforms
        self.transforms = transforms.Compose([
            transforms.Resize((50, 50)),
            transforms.ToTensor(),
            transforms.Normalize([0.485, 0.456, 0.406], [0.229, 0.224, 0.225])
        ])

    def predict(self, image: Image.Image):
        # Preprocess the image
        image_tensor = self.transforms(image).unsqueeze(0).to(self.device)

        # Perform inference
        with torch.no_grad():
            outputs = self.model(image_tensor)
            probabilities = torch.nn.functional.softmax(outputs, dim=1).cpu().numpy().tolist()
            label = int(torch.argmax(outputs, dim=1))
        return label, probabilities
