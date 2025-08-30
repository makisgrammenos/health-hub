#!/usr/bin/env python3
"""
Inference for BI-RADS A–D using MONAI TorchVisionFCModel (inception_v3).

Features:
- Uses the same model wrapper as the MONAI bundle → perfect checkpoint compatibility.
- Prints predictions directly to the terminal.
- No CSV, no output files, no saving.

Usage:
  python inference_monai_bundle.py --inputs sample_data/A/sample_A1.jpg --checkpoint ./models/model.pt
  python inference_monai_bundle.py --inputs "./sample_data/*/*.jpg"
"""

import argparse
import json
import os
import sys
import glob
from pathlib import Path
from typing import List, Dict, Any

import torch
from monai.data import Dataset, DataLoader
from monai.transforms import (
    Compose, LoadImaged, EnsureChannelFirstd, ScaleIntensityd, Resized, Activationsd
)
from monai.networks.nets.torchvision_fc import TorchVisionFCModel

IMG_EXTS = {".jpg", ".jpeg", ".png", ".bmp", ".tif", ".tiff"}
BI_RADS_MAP = {0: "A", 1: "B", 2: "C", 3: "D"}


# ---------------------------------------------------
# Resolve Inputs
# ---------------------------------------------------
def resolve_inputs(inputs_arg: str, split_name: str | None) -> List[str]:
    p = Path(inputs_arg)
    if p.is_file():
        suf = p.suffix.lower()
        if suf == ".json":
            with open(p, "r", encoding="utf-8") as f:
                payload = json.load(f)
            if split_name is None:
                for k in ("Test", "test", "Val", "Validation", "infer", "Inference"):
                    if k in payload:
                        split_name = k
                        break
            if split_name is None or split_name not in payload:
                raise ValueError(f"Split '{split_name}' not found in {p}. Available: {list(payload.keys())}")
            items = payload[split_name]
            if not items:
                return []
            if isinstance(items[0], str):
                return items
            if isinstance(items[0], dict):
                if "image" in items[0]:
                    return [d["image"] for d in items if "image" in d]
                raise ValueError("JSON items are dicts but missing 'image' key.")
            raise ValueError("Unsupported JSON format.")
        if suf in IMG_EXTS:
            return [str(p)]
        with open(p, "r", encoding="utf-8") as f:
            return [ln.strip() for ln in f if ln.strip()]

    if p.is_dir():
        paths: List[str] = []
        for ext in IMG_EXTS:
            paths.extend(sorted(str(q) for q in p.rglob(f"*{ext}")))
        return paths

    matches = glob.glob(inputs_arg)
    if matches:
        return sorted(matches)

    raise ValueError(f"Could not resolve inputs from: {inputs_arg}")


# ---------------------------------------------------
# Transforms
# ---------------------------------------------------
def build_transforms() -> Compose:
    return Compose([
        LoadImaged(keys="image"),
        EnsureChannelFirstd(keys="image"),
        ScaleIntensityd(keys="image", minv=0.0, maxv=1.0),
        Resized(keys="image", spatial_size=(299, 299)),
    ])


def make_dataset(image_paths: List[str], transforms: Compose) -> Dataset:
    data = [{"image": p} for p in image_paths]
    return Dataset(data=data, transform=transforms)


# ---------------------------------------------------
# Load Model
# ---------------------------------------------------
def load_model(device: torch.device, num_classes: int, checkpoint_path: str | None) -> torch.nn.Module:
    model = TorchVisionFCModel(
        model_name="inception_v3",
        num_classes=num_classes,
        pool=None,
        use_conv=False,
        bias=True,
        pretrained=False,
    ).to(device)
    model.eval()

    if checkpoint_path and os.path.isfile(checkpoint_path):
        ckpt = torch.load(checkpoint_path, map_location="cpu")
        state = ckpt.get("model", ckpt)
        missing, unexpected = model.load_state_dict(state, strict=False)
        if missing or unexpected:
            print(f"[warn] Loaded with non-strict matching.\nMissing: {missing}\nUnexpected: {unexpected}", file=sys.stderr)
        else:
            print("[info] Checkpoint loaded cleanly.")
    else:
        print(f"[info] No checkpoint at {checkpoint_path}; using randomly initialized weights.", file=sys.stderr)

    return model


# ---------------------------------------------------
# Inference
# ---------------------------------------------------
def infer(
    model: torch.nn.Module,
    loader: DataLoader,
    device: torch.device,
) -> List[Dict[str, Any]]:
    post = Activationsd(keys="pred", sigmoid=True)
    results: List[Dict[str, Any]] = []
    amp = device.type == "cuda"

    with torch.no_grad():
        for batch in loader:
            imgs = batch["image"].to(device=device, dtype=torch.float32)
            paths = batch.get("image_meta_dict", {}).get("filename_or_obj", None)
            if paths is None:
                paths = batch.get("path", None)

            with torch.cuda.amp.autocast(enabled=amp):
                logits = model(imgs)
                pred = {"pred": logits}

            pred = post(pred)
            probs = pred["pred"].float().cpu().numpy()  # (B,4)

            for i in range(probs.shape[0]):
                p = probs[i].tolist()
                pred_idx = int(max(range(4), key=lambda j: p[j]))
                results.append({
                    "image_path": paths[i] if isinstance(paths, list) else None,
                    "prob_A": p[0], "prob_B": p[1], "prob_C": p[2], "prob_D": p[3],
                    "pred_idx": pred_idx,
                    "pred_label": BI_RADS_MAP.get(pred_idx, str(pred_idx)),
                })

    return results


# ---------------------------------------------------
# Main
# ---------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="MONAI-bundle-accurate inference with TorchVisionFCModel (inception_v3)")
    ap.add_argument("--inputs", required=True, help="Image path | directory | glob | JSON dataset | text file of paths")
    ap.add_argument("--split", default="Test", help="Split name if --inputs is JSON (default: Test)")
    ap.add_argument("--checkpoint", default="./models/model.pt", help="Path to model checkpoint")
    ap.add_argument("--batch-size", type=int, default=4)
    ap.add_argument("--num-workers", type=int, default=2)
    args = ap.parse_args()

    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    print(f"[info] Using device: {device}")

    try:
        paths = resolve_inputs(args.inputs, args.split)
    except Exception as e:
        print(f"[error] {e}", file=sys.stderr)
        sys.exit(1)

    if not paths:
        print("[error] No images resolved from inputs.", file=sys.stderr)
        sys.exit(1)

    print(f"[info] Found {len(paths)} image(s).")

    transforms = build_transforms()
    ds = make_dataset(paths, transforms)
    loader = DataLoader(
        ds,
        batch_size=args.batch_size,
        shuffle=False,
        num_workers=args.num_workers,
        pin_memory=(device.type == "cuda"),
    )

    model = load_model(device, num_classes=4, checkpoint_path=args.checkpoint)

    results = infer(model, loader, device)

    # Print results to terminal only
    for r in results:
        print(f"\nImage: {r['image_path']}")
        print(f"  Probabilities: A={r['prob_A']:.4f}, B={r['prob_B']:.4f}, C={r['prob_C']:.4f}, D={r['prob_D']:.4f}")
        print(f"  Prediction: {r['pred_label']} (class {r['pred_idx']})")


if __name__ == "__main__":
    main()
