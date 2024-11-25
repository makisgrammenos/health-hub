import torch
from monai.transforms import Activations, AsDiscrete

def preprocess_data(data_dict, transforms):
    """
    Apply preprocessing transforms to the input data dictionary.
    """
    data = transforms(data_dict)
    inputs = data["image"].unsqueeze(0).to(torch.device("cuda:0" if torch.cuda.is_available() else "cpu"))
    return inputs

def postprocess_predictions(predictions):
    """
    Postprocess the model predictions to generate segmentation labels.
    """
    activations = Activations(sigmoid=True)
    discrete = AsDiscrete(threshold=0.5)
    predictions = activations(predictions)
    predictions = discrete(predictions)

    # Convert predictions to numpy
    pred_image = predictions[0].cpu().numpy()

    return {"segmentation": pred_image.tolist()}
