# models/pancreas_seg/model.py
import os
import io
import uuid
import asyncio
import tempfile
from typing import List, Tuple, Dict

import numpy as np
import nibabel as nib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

import torch
# match your script: stick to GPU 0 if present
if torch.cuda.is_available():
    torch.cuda.set_device(0)

# ---- PyTorch 2.6 checkpoint compatibility (trusted checkpoints only) ----
_real_torch_load = torch.load
def _torch_load_unsafe(*args, **kwargs):
    kwargs.setdefault("weights_only", False)
    return _real_torch_load(*args, **kwargs)
torch.load = _torch_load_unsafe
try:
    from numpy.core.multiarray import _reconstruct as _np_reconstruct
    torch.serialization.add_safe_globals([_np_reconstruct])
except Exception:
    pass
# -------------------------------------------------------------------------

from monai.bundle import ConfigParser
from monai.transforms import (
    Compose, LoadImaged, EnsureChannelFirstd, Orientationd,
    Spacingd, ScaleIntensityRanged, EnsureTyped
)
from monai.inferers import SlidingWindowInferer

# Concurrency control for GPU
GPU_LOCK = asyncio.Semaphore(int(os.getenv("MAX_CONCURRENT_GPU", "1")))

# Robustly locate inference.yaml next to this file by default
BUNDLE_CONFIG = os.getenv(
    "PANCREAS_BUNDLE_CONFIG",
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "inference.yaml"),
)

# Simple overlay palette: 0 bg (transparent), 1 pancreas (green), 2 lesion (red)
SEG_COLORS = [
    (0, 0, 0, 0),
    (0, 1, 0, 0.6),
    (1, 0, 0, 0.6),
]
SEG_CMAP = ListedColormap(SEG_COLORS, name="pancreas_cmap")
SEG_NORM = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], ncolors=len(SEG_COLORS))


def _overlay_png(mri_slice: np.ndarray, seg_slice: np.ndarray) -> str:
    """Return a data-url PNG of CT + segmentation overlay (uses basic WL/WW-ish normalization)."""
    # Normalize CT slice to [0,1] for display (robust percentiles)
    if np.isfinite(mri_slice).any():
        vmin, vmax = np.percentile(mri_slice, (1, 99))
        mri_norm = np.clip((mri_slice - vmin) / (vmax - vmin + 1e-8), 0, 1)
    else:
        mri_norm = np.zeros_like(mri_slice, dtype=np.float32)

    fig = plt.figure(figsize=(8, 8), dpi=100)
    ax = plt.axes([0, 0, 1, 1])
    ax.imshow(mri_norm, cmap="gray", interpolation="nearest")
    ax.imshow(seg_slice, cmap=SEG_CMAP, norm=SEG_NORM, alpha=0.7, interpolation="nearest")
    ax.axis("off")
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=100)
    plt.close(fig)
    buf.seek(0)
    import base64
    return "data:image/png;base64," + base64.b64encode(buf.read()).decode("utf-8")


async def predict_single_nifti(input_path: str) -> Tuple[str, List[Dict], Dict[str, int]]:
    """
    Run single-image segmentation using the same steps as your working script.
    Returns:
      - mask_path: saved NIfTI mask path (prediction.nii.gz, affine from original image)
      - overlays: [{slice_index, image(data-url)} ...] built from resampled CT + mask
      - volumes: voxel counts
    """
    # --- 1) Load bundle + network (exactly like your script) ---
    parser = ConfigParser()
    parser.read_config(BUNDLE_CONFIG)
    # ensure bundle_root resolves YAML paths like "$@bundle_root + '/models/model.pt'"
    parser.config["bundle_root"] = os.path.dirname(BUNDLE_CONFIG)
    parser.parse()

    device = parser.get_parsed_content("device")
    network = parser.get_parsed_content("network").to(device).eval()

    # --- 2) Load checkpoint (using bundle's CheckpointLoader.load_path) ---
    def _load_weights_from_bundle():
        if "checkpointloader" not in parser.config:
            print("[pancreas] No checkpointloader in config; using randomly initialized weights.")
            return
        cl = parser.get_parsed_content("checkpointloader")
        load_path = getattr(cl, "load_path", None)
        if not load_path or not os.path.isfile(load_path):
            print(f"[pancreas][WARN] Checkpoint path not found: {load_path}. Skipping weights.")
            return
        print(f"[pancreas] Loading weights from: {load_path}")
        state = torch.load(load_path, map_location=device)
        loaded = False
        for key in ("state_dict", "net", "model", "model_state", "network"):
            if isinstance(state, dict) and key in state and isinstance(state[key], dict):
                network.load_state_dict(state[key], strict=False)
                loaded = True
                break
        if not loaded and isinstance(state, dict):
            try:
                network.load_state_dict(state, strict=False)
                loaded = True
            except Exception:
                pass
        print("[pancreas] Weights loaded." if loaded else "[pancreas][WARN] Unrecognized checkpoint format.")
    _load_weights_from_bundle()

    # --- 3) Preprocessing (mirrors your script) ---
    preprocess = Compose([
        LoadImaged(keys="image"),
        EnsureChannelFirstd(keys="image"),
        Orientationd(keys="image", axcodes="RAS"),
        Spacingd(keys="image", pixdim=(1, 1, 1), mode="bilinear"),
        ScaleIntensityRanged(keys="image", a_min=-87, a_max=199, b_min=0, b_max=1, clip=True),
        EnsureTyped(keys="image", track_meta=True),
    ])

    # --- 4) Infer (same ROI/overlap as your script) ---
    inferer = SlidingWindowInferer(roi_size=(96, 96, 96), sw_batch_size=4, overlap=0.625)

    out_dir = tempfile.gettempdir()
    base_noext = os.path.basename(input_path).replace(".nii.gz", "").replace(".nii", "")
    mask_path = os.path.join(out_dir, "prediction.nii.gz")  # same filename as your script
    resampled_img_path = os.path.join(out_dir, "image_resampled.nii.gz")

    async with GPU_LOCK:
        data = preprocess({"image": input_path})
        img = data["image"].unsqueeze(0).to(device)  # [B,C,D,H,W]

        with torch.no_grad():
            logits = inferer(inputs=img, network=network)
            softmax = torch.softmax(logits, dim=1)
            pred = torch.argmax(softmax, dim=1).cpu().numpy()[0]  # [D,H,W], uint8 later

        # Save resampled image (to match mask dims visually)
        img_mt = data["image"]
        resampled_vol = img_mt.cpu().numpy()[0]
        # Affine on resampled grid (robust fallbacks)
        try:
            affine_resampled = img_mt.affine.numpy()
        except Exception:
            try:
                affine_resampled = data["image_meta_dict"]["affine"].numpy()
            except Exception:
                spacing = (getattr(img_mt, "meta", {}) or {}).get("spacing", None)
                if spacing is not None:
                    sx, sy, sz = map(float, spacing)
                    affine_resampled = np.diag([sx, sy, sz, 1.0]).astype(np.float32)
                else:
                    affine_resampled = np.eye(4, dtype=np.float32)
        nib.save(nib.Nifti1Image(resampled_vol.astype(np.float32), affine_resampled), resampled_img_path)

        # Save mask (your script uses original affine)
        ref_img = nib.load(input_path)
        nib.save(nib.Nifti1Image(pred.astype(np.uint8), affine=ref_img.affine), mask_path)

    # --- 5) Build overlays from resampled CT + mask (dims match) ---
    overlays = []
    D = int(pred.shape[0])
    for z in range(D):
        overlays.append({
            "slice_index": z,
            "image": _overlay_png(resampled_vol[z, :, :], pred[z, :, :]),
        })

    # --- 6) Simple voxel counts ---
    volumes = {
        "pancreas_voxels": int((pred == 1).sum()),
        "lesion_voxels": int((pred == 2).sum()),
        "total_seg_voxels": int((pred > 0).sum()),
    }

    return mask_path, overlays, volumes
