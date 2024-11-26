import os
import torch
from monai.networks.nets import SegResNet
from monai.transforms import Compose, LoadImaged, EnsureChannelFirstd, NormalizeIntensityd, ConcatItemsd
from monai.inferers import SlidingWindowInferer
from .config import MODEL_PATH, DEVICE
from .utils import preprocess_data, postprocess_predictions

# Load the model
def load_model():
    """
    Load the SegResNet model for brain tumor segmentation.
    """
    if not os.path.exists(MODEL_PATH):
        raise RuntimeError(f"Model checkpoint not found: {MODEL_PATH}")

    model = SegResNet(
        blocks_down=[1, 2, 2, 4],
        blocks_up=[1, 1, 1],
        init_filters=16,
        in_channels=4,
        out_channels=3,
        dropout_prob=0.2,
    ).to(DEVICE)

    checkpoint = torch.load(MODEL_PATH, map_location=DEVICE)
    model.load_state_dict(checkpoint)
    model.eval()
    return model

# Perform segmentation
def segment_brain_tumor(file_paths, model):
    """
    Perform segmentation on the provided MRI modalities.
    """
    # Preprocessing pipeline
    transforms = Compose([
        LoadImaged(keys=list(file_paths.keys())),
        EnsureChannelFirstd(keys=list(file_paths.keys())),
        ConcatItemsd(keys=list(file_paths.keys()), name="image", dim=0),
        NormalizeIntensityd(keys="image", nonzero=True, channel_wise=True),
    ])

    # Prepare data dictionary
    data_dict = {key: file_paths[key] for key in file_paths}

    # Preprocess the data
    inputs = preprocess_data(data_dict, transforms)

    # Sliding window inference
    inferer = SlidingWindowInferer(roi_size=(240, 240, 160), sw_batch_size=1, overlap=0.5)
    with torch.no_grad():
        predictions = inferer(inputs, model)

    # Postprocess predictions
    return postprocess_predictions(predictions)
