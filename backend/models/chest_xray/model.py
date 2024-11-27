import torch
import torchxrayvision as xrv

class ChestXRayModel:
    def __init__(self):
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.model = xrv.models.DenseNet(weights="densenet121-res224-all")
        self.model = self.model.to(self.device)
        self.model.eval()

    def predict(self, img_tensor):
        with torch.no_grad():
            return self.model(img_tensor)

    def get_pathologies(self):
        return self.model.pathologies
