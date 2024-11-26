import torch

# Model configuration
MODEL_CHECKPOINT_PATH = "models/brain_tumor_segmentation/model.pt"
DEVICE = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
