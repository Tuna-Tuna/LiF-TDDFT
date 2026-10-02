"""In-domain interpolation of calculated detachment probabilities only."""

from __future__ import annotations

import numpy as np


class OutOfDomainError(ValueError):
    pass


class DetachmentProbabilityInterpolator:
    def __init__(self, heights: np.ndarray, velocities: np.ndarray, probabilities: np.ndarray):
        self.heights = np.asarray(heights, dtype=float)
        self.velocities = np.asarray(velocities, dtype=float)
        self.values = np.asarray(probabilities, dtype=float)
        for axis in (self.heights, self.velocities):
            if axis.ndim != 1 or axis.size < 2 or np.any(~np.isfinite(axis)):
                raise ValueError("probability grid axes need at least two finite nodes")
            if np.any(np.diff(axis) <= 0):
                raise ValueError("probability grid axes must be strictly increasing")
        if self.values.shape != (self.heights.size, self.velocities.size):
            raise ValueError("probabilities shape must be (n_heights, n_velocities)")
        if np.any(~np.isfinite(self.values)):
            raise ValueError("missing calculated probability nodes are not filled")
        if np.any(self.values < 0.0) or np.any(self.values > 1.0):
            raise ValueError("probabilities must lie in [0,1]; clipping or rescaling is forbidden")

    def __call__(self, height: float, velocity: float) -> float:
        h, v = float(height), float(velocity)
        if not (self.heights[0] <= h <= self.heights[-1]
                and self.velocities[0] <= v <= self.velocities[-1]):
            raise OutOfDomainError(f"event (h={h}, v={v}) lies outside the TD probability database")
        ih, iv = np.flatnonzero(self.heights == h), np.flatnonzero(self.velocities == v)
        if ih.size and iv.size:
            return float(self.values[ih[0], iv[0]])
        i = min(max(int(np.searchsorted(self.heights, h) - 1), 0), self.heights.size - 2)
        j = min(max(int(np.searchsorted(self.velocities, v) - 1), 0), self.velocities.size - 2)
        h0, h1 = self.heights[i:i + 2]
        v0, v1 = self.velocities[j:j + 2]
        th, tv = (h - h0) / (h1 - h0), (v - v0) / (v1 - v0)
        q00, q01 = self.values[i, j], self.values[i, j + 1]
        q10, q11 = self.values[i + 1, j], self.values[i + 1, j + 1]
        value = (1 - th) * ((1 - tv) * q00 + tv * q01) + th * ((1 - tv) * q10 + tv * q11)
        if not 0.0 <= value <= 1.0:
            raise ValueError("interpolated probability outside [0,1]; clipping is forbidden")
        return float(value)
