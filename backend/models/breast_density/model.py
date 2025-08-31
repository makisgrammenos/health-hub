import os
from typing import Tuple, List, Dict, Any

import numpy as np
import torch
from PIL import Image

from monai.transforms import Compose, EnsureChannelFirstd, ScaleIntensityd, Resized
from monai.networks.nets.torchvision_fc import TorchVisionFCModel

# Class order must match training / bundle / CLI script
CLASS_ORDER: List[str] = ["A", "B", "C", "D"]  # indices 0..3


class BiradsDensityModel:
    def __init__(self, checkpoint_path: str, device: str | None = None, debug: bool = False):
        self.device = torch.device(device or ("cuda:0" if torch.cuda.is_available() else "cpu"))
        self.debug = debug

        # Match the CLI script: ImageNet init, eval() before loading checkpoint
        self.model = TorchVisionFCModel(
            model_name="inception_v3",
            num_classes=len(CLASS_ORDER),
            pool=None,
            use_conv=False,
            bias=True,
            pretrained=True,   # <- critical (same as script)
        ).to(self.device)
        self.model.eval()      # <- set eval BEFORE loading weights

        if not os.path.isfile(checkpoint_path):
            raise FileNotFoundError(f"Checkpoint not found: {checkpoint_path}")

        if self.device.type == "cuda":
            torch.backends.cudnn.benchmark = True  # parity with script

        # Load checkpoint mapped to target device (parity with script)
        print(f"[birads] loading checkpoint: {checkpoint_path}")
        ckpt = torch.load(checkpoint_path, map_location=self.device)

        # Handle common checkpoint dict layouts
        if isinstance(ckpt, dict):
            if "model" in ckpt:
                state = ckpt["model"]
            elif "state_dict" in ckpt:
                state = ckpt["state_dict"]
            elif "network" in ckpt:
                state = ckpt["network"]
            else:
                state = ckpt
        else:
            state = ckpt

        # Strict first, then fallback non-strict (exactly like the script)
        try:
            self.model.load_state_dict(state, strict=True)
            print("[birads] checkpoint loaded (strict)")
        except Exception as e:
            print(f"[birads] strict load failed, trying non-strict: {e}")
            missing, unexpected = self.model.load_state_dict(state, strict=False)
            if missing:
                print(f"[birads] missing keys: {missing}")
            if unexpected:
                print(f"[birads] unexpected keys: {unexpected}")

        # Dict-style transforms with channel_dim=2 (HWC -> CHW), exactly like script
        self.tx = Compose([
            EnsureChannelFirstd(keys="image", channel_dim=2),
            ScaleIntensityd(keys="image", minv=0.0, maxv=1.0),
            Resized(keys="image", spatial_size=(299, 299)),
        ])

    @torch.inference_mode()
    def predict(self, pil_image: Image.Image) -> Tuple[int, List[float]]:
        """
        Returns (pred_idx, probs[4]) with sigmoid activation (bundle/script behavior).
        """
        # Normalize PIL mode; keep L or RGB. Other modes -> RGB for safety.
        if pil_image.mode not in ("RGB", "L"):
            pil_image = pil_image.convert("RGB")

        arr = np.array(pil_image)  # HxW or HxWx3 (uint8)

        # Dict-style pipeline as in the script (we already have the pixel array; no LoadImaged here)
        sample: Dict[str, Any] = {"image": arr}
        sample = self.tx(sample)     # -> {'image': CHW float32 in [0,1]}
        x = sample["image"]

        # Enforce 3 channels per metadata (grayscale -> 3ch copy), same as script intent
        if isinstance(x, np.ndarray):
            if x.ndim == 2:
                x = x[None, ...]
            if x.shape[0] == 1:
                x = np.repeat(x, 3, axis=0)
            elif x.shape[0] > 3:
                x = x[:3]
            x = torch.from_numpy(x).float()
        else:
            # Tensor path
            if x.ndim == 2:
                x = x.unsqueeze(0)
            if x.shape[0] == 1:
                x = x.repeat(3, 1, 1)
            elif x.shape[0] > 3:
                x = x[:3]
            x = x.float()

        x = x.unsqueeze(0).to(self.device)  # 1x3x299x299

        if self.debug:
            print(f"[birads] input tensor: {tuple(x.shape)} "
                  f"range[{x.min().item():.3f},{x.max().item():.3f}] device={x.device}")

        amp = self.device.type == "cuda"
        with torch.cuda.amp.autocast(enabled=amp):
            out = self.model(x)
            # Handle potential tuple (aux outputs); script expects main logits
            logits = out[0] if isinstance(out, tuple) else out

        if logits.shape[-1] != len(CLASS_ORDER):
            raise RuntimeError(f"Unexpected logits shape {tuple(logits.shape)}; expected last dim {len(CLASS_ORDER)}")

        probs_t = torch.sigmoid(logits)  # script uses sigmoid via Activationsd
        probs = probs_t.squeeze(0).tolist()
        pred_idx = int(np.argmax(probs))  # same as script

        return pred_idx, probs

    @staticmethod
    def format_response(filename: str, pred_idx: int, probs: List[float]) -> Dict[str, Any]:
        return {
            "filename": filename,
            "prediction_index": pred_idx,
            "prediction_label": CLASS_ORDER[pred_idx],
            "probabilities": {
                "A": float(probs[0]),
                "B": float(probs[1]),
                "C": float(probs[2]),
                "D": float(probs[3]),
            },
        }
