import torch
from monai.networks.nets import SegResNet
from monai.transforms import Compose, Activations, AsDiscrete
from monai.inferers import SlidingWindowInferer
import numpy as np


class BrainTumorSegmentationModel:
    def __init__(self, model_path: str, device: torch.device):
        self.device = device
        self.model = SegResNet(
            blocks_down=[1, 2, 2, 4],
            blocks_up=[1, 1, 1],
            init_filters=16,
            in_channels=4,
            out_channels=3,
            dropout_prob=0.2
        ).to(self.device)

        self.model.load_state_dict(torch.load(model_path, map_location=self.device))
        self.model.eval()

        self.inferer = SlidingWindowInferer(
            roi_size=(240, 240, 160),
            sw_batch_size=1,
            overlap=0.5
        )

        self.post_process = Compose([
            Activations(sigmoid=True),
            AsDiscrete(threshold=0.5)
        ])

    def predict(self, input_tensor: torch.Tensor) -> np.ndarray:
        with torch.no_grad():
            prediction = self.inferer(input_tensor.to(self.device), self.model)
            prediction = self.post_process(prediction)
        return prediction[0].cpu().numpy()
