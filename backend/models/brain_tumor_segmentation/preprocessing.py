from monai.transforms import (
    Compose,
    LoadImaged,
    EnsureChannelFirstd,
    NormalizeIntensityd,
    ConcatItemsd
)


def get_preprocessing_transforms(modalities):
    return Compose([
        LoadImaged(keys=modalities),
        EnsureChannelFirstd(keys=modalities),
        ConcatItemsd(keys=modalities, name="image", dim=0),
        NormalizeIntensityd(keys="image", nonzero=True, channel_wise=True)
    ])
