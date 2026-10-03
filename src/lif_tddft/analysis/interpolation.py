"""PRA height-only shape-preserving interpolation of Pdet, without extrapolation."""
import numpy as np
from scipy.interpolate import PchipInterpolator


class OutOfDomainError(ValueError):
    pass


class DetachmentProbabilityInterpolator:
    def __init__(self, heights, velocities, probabilities):
        self.heights = np.asarray(heights, dtype=float)
        self.velocities = np.asarray(velocities, dtype=float)
        self.values = np.asarray(probabilities, dtype=float)
        for axis, minimum in ((self.heights, 2), (self.velocities, 1)):
            if axis.ndim != 1 or axis.size < minimum or np.any(~np.isfinite(axis)):
                raise ValueError('probability axes require finite recorded nodes')
            if np.any(np.diff(axis) <= 0):
                raise ValueError('probability grid axes must be strictly increasing')
        if self.values.shape != (self.heights.size, self.velocities.size):
            raise ValueError('probabilities shape must be (n_heights, n_velocities)')
        if np.any(~np.isfinite(self.values)):
            raise ValueError('missing calculated probability nodes are not filled')
        if np.any(self.values < 0) or np.any(self.values > 1):
            raise ValueError('probabilities must lie in [0,1]; clipping or rescaling is forbidden')
        self.curves = PchipInterpolator(self.heights, self.values, axis=0, extrapolate=False)

    def __call__(self, height, velocity):
        h, v = float(height), float(velocity)
        if not np.isfinite(h) or not self.heights[0] <= h <= self.heights[-1]:
            raise OutOfDomainError(f'height {h} lies outside the calculated detachment grid')
        matches = np.flatnonzero(np.isclose(self.velocities, v, rtol=0, atol=1e-12))
        if matches.size != 1:
            raise OutOfDomainError('velocity must match one calculated velocity; no velocity interpolation')
        j = int(matches[0])
        ih = np.flatnonzero(self.heights == h)
        if ih.size:
            return float(self.values[ih[0], j])
        value = float(self.curves(h)[j])
        if not np.isfinite(value) or not 0 <= value <= 1:
            raise ValueError('interpolated probability outside [0,1]; clipping is forbidden')
        return value
