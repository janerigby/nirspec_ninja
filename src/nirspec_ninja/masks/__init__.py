"""Custom NIRSpec masks for NSClean/CleanFlickerNoise."""
from .custom_nirspec_mask_2 import make_custom_NIRSpec_mask, refine_mask_from_rate, refine_mask_file_from_rate

__all__ = ["make_custom_NIRSpec_mask", "refine_mask_from_rate", "refine_mask_file_from_rate"]
