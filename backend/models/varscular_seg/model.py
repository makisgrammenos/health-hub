import torch
import numpy as np
from monai.networks.nets import UNet
from monai.transforms import Compose, Activations, AsDiscrete
from monai.inferers import sliding_window_inference
import logging

logger = logging.getLogger(__name__)


class MonaiUNet2DSegmentationModel:
    """
    MONAI UNet 2D Segmentation Model for DICOM processing
    """
    
    def __init__(self, model_path: str, device: torch.device):
        """
        Initialize the UNet 2D model
        
        Args:
            model_path: Path to the model checkpoint
            device: torch device (cuda/cpu)
        """
        self.device = device
        
        # Initialize UNet architecture (matching your config)
        self.model = UNet(
            spatial_dims=2,
            in_channels=1,
            out_channels=4,  # 4 classes for segmentation
            channels=[16, 32, 64, 128, 256],
            strides=[2, 2, 2, 2],
            num_res_units=2
        ).to(self.device)
        
        # Load checkpoint
        logger.info(f"Loading model from {model_path}")
        checkpoint = torch.load(model_path, map_location=self.device)
        
        # Handle different checkpoint formats
        if 'model' in checkpoint:
            self.model.load_state_dict(checkpoint['model'])
        elif 'state_dict' in checkpoint:
            self.model.load_state_dict(checkpoint['state_dict'])
        else:
            self.model.load_state_dict(checkpoint)
        
        self.model.eval()
        logger.info(f"Model loaded successfully on {device}")
        
        # Post-processing transforms
        self.post_process = Compose([
            Activations(softmax=True),
            AsDiscrete(argmax=True)
        ])
        
        # Sliding window parameters
        self.roi_size = (256, 256)
        self.sw_batch_size = 1
        self.overlap = 0.25
    
    def predict(self, input_tensor: torch.Tensor) -> np.ndarray:
        """
        Run inference on input tensor
        
        Args:
            input_tensor: Input tensor with shape (B, C, H, W) for 2D
                         or (B, C, D, H, W) for 3D volumes
        
        Returns:
            Segmentation prediction as numpy array
        """
        with torch.no_grad():
            # Check input dimensions
            if len(input_tensor.shape) == 4:  # 2D: (batch, channel, height, width)
                prediction = sliding_window_inference(
                    inputs=input_tensor.to(self.device),
                    roi_size=self.roi_size,
                    sw_batch_size=self.sw_batch_size,
                    predictor=self.model,
                    overlap=self.overlap
                )
            elif len(input_tensor.shape) == 5:  # 3D: (batch, channel, depth, height, width)
                # Process slice by slice for 3D volumes
                b, c, d, h, w = input_tensor.shape
                predictions = []
                
                for slice_idx in range(d):
                    slice_input = input_tensor[:, :, slice_idx:slice_idx+1, :, :]
                    # Reshape to 2D for processing
                    slice_input = slice_input.squeeze(2).unsqueeze(0)
                    
                    slice_pred = sliding_window_inference(
                        inputs=slice_input.to(self.device),
                        roi_size=self.roi_size,
                        sw_batch_size=self.sw_batch_size,
                        predictor=self.model,
                        overlap=self.overlap
                    )
                    predictions.append(slice_pred)
                
                # Stack predictions back to 3D
                prediction = torch.stack(predictions, dim=2)
            else:
                raise ValueError(f"Unexpected input shape: {input_tensor.shape}")
            
            # Apply post-processing
            prediction = self.post_process(prediction)
            
        return prediction[0].cpu().numpy()  # Remove batch dimension and convert to numpy
    
    def predict_single_slice(self, input_slice: torch.Tensor) -> np.ndarray:
        """
        Run inference on a single 2D slice
        
        Args:
            input_slice: Input tensor with shape (C, H, W)
        
        Returns:
            Segmentation prediction as numpy array
        """
        # Add batch dimension
        input_tensor = input_slice.unsqueeze(0).to(self.device)
        
        with torch.no_grad():
            prediction = sliding_window_inference(
                inputs=input_tensor,
                roi_size=self.roi_size,
                sw_batch_size=self.sw_batch_size,
                predictor=self.model,
                overlap=self.overlap
            )
            
            # Apply post-processing
            prediction = self.post_process(prediction)
        
        return prediction[0].cpu().numpy()  # Remove batch dimension