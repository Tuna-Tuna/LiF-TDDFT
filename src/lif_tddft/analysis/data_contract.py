"""Hard gates preventing probability-data scaling or silent correction."""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np


def assert_no_scaling_policy(metadata: Mapping[str, object]) -> None:
    """Require an explicit identity-transform policy for figure/yield data."""

    required = {
        "normalization": "none",
        "scale_factor": 1.0,
        "value_transform": "none",
        "velocity_specific_adjustments": "none",
    }
    missing = [key for key in required if key not in metadata]
    if missing:
        raise ValueError(f"no-scaling metadata is incomplete: {missing}")
    if str(metadata["normalization"]).strip().lower() != "none":
        raise ValueError("probability normalization is forbidden")
    if float(metadata["scale_factor"]) != 1.0:
        raise ValueError("probability scale_factor must be exactly 1.0")
    if str(metadata["value_transform"]).strip().lower() != "none":
        raise ValueError("probability value transforms are forbidden")
    if str(metadata["velocity_specific_adjustments"]).strip().lower() != "none":
        raise ValueError("velocity-specific probability adjustments are forbidden")


def validate_probability_values(values: np.ndarray, *, label: str) -> np.ndarray:
    """Validate probabilities without clipping, normalization, or rescaling."""

    array = np.asarray(values, dtype=float)
    if np.any(~np.isfinite(array)):
        raise ValueError(f"{label} contains non-finite values")
    if np.any(array < 0.0) or np.any(array > 1.0):
        raise ValueError(f"{label} lies outside [0,1]; clipping is forbidden")
    return array


def assert_values_identical(source: np.ndarray, exported: np.ndarray, *, label: str) -> None:
    """Require bitwise-identical numeric arrays at a data-export boundary."""

    source_array = np.asarray(source, dtype=float)
    exported_array = np.asarray(exported, dtype=float)
    if source_array.shape != exported_array.shape or not np.array_equal(
        source_array, exported_array
    ):
        raise ValueError(f"{label} changed across export; scaling or rewriting is forbidden")
