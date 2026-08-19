"""Grid interpolation with explicit domain and value contracts."""

from __future__ import annotations

import numpy as np


class OutOfDomainError(ValueError):
    pass


class BoundedGridInterpolator:
    def __init__(
        self,
        heights: np.ndarray,
        velocities: np.ndarray,
        values: np.ndarray,
        *,
        value_bounds: tuple[float, float] | None = (0.0, 1.0),
    ):
        self.heights = np.asarray(heights, dtype=float)
        self.velocities = np.asarray(velocities, dtype=float)
        self.values = np.asarray(values, dtype=float)
        if self.values.shape != (self.heights.size, self.velocities.size):
            raise ValueError("values shape must be (n_heights, n_velocities)")
        if np.any(np.diff(self.heights) <= 0) or np.any(np.diff(self.velocities) <= 0):
            raise ValueError("grid axes must be strictly increasing")
        if np.any(~np.isfinite(self.values)):
            raise ValueError("missing/non-finite production data cannot be interpolated")
        if value_bounds is not None:
            lower, upper = map(float, value_bounds)
            if not lower < upper:
                raise ValueError("value_bounds must be strictly increasing")
            if np.any(self.values < lower) or np.any(self.values > upper):
                raise ValueError(
                    "input values lie outside value_bounds; clipping or rescaling is forbidden"
                )
            self.value_bounds: tuple[float, float] | None = (lower, upper)
        else:
            self.value_bounds = None

    def __call__(self, height: float, velocity: float) -> float:
        h = float(height)
        v = float(velocity)
        inside = self.heights[0] <= h <= self.heights[-1] and self.velocities[0] <= v <= self.velocities[-1]
        if not inside:
            raise OutOfDomainError(f"event (h={h}, v={v}) lies outside the TD database")
        i = min(max(int(np.searchsorted(self.heights, h) - 1), 0), self.heights.size - 2)
        j = min(max(int(np.searchsorted(self.velocities, v) - 1), 0), self.velocities.size - 2)
        h0, h1 = self.heights[i : i + 2]
        v0, v1 = self.velocities[j : j + 2]
        th = (h - h0) / (h1 - h0)
        tv = (v - v0) / (v1 - v0)
        q00, q01 = self.values[i, j], self.values[i, j + 1]
        q10, q11 = self.values[i + 1, j], self.values[i + 1, j + 1]
        value = (1 - th) * ((1 - tv) * q00 + tv * q01) + th * ((1 - tv) * q10 + tv * q11)
        if self.value_bounds is not None:
            lower, upper = self.value_bounds
            if value < lower or value > upper:
                raise ValueError(
                    "interpolated value lies outside value_bounds; clipping is forbidden"
                )
        return float(value)

    def leave_one_out_error(self) -> float:
        """Inverse-distance leave-one-grid-point-out RMS error."""
        hh, vv = np.meshgrid(self.heights, self.velocities, indexing="ij")
        points = np.column_stack((hh.ravel(), vv.ravel()))
        values = self.values.ravel()
        h_scale = max(float(np.ptp(self.heights)), 1.0e-15)
        v_scale = max(float(np.ptp(self.velocities)), 1.0e-15)
        errors = []
        for index, point in enumerate(points):
            mask = np.arange(points.shape[0]) != index
            delta = points[mask] - point
            distance2 = (delta[:, 0] / h_scale) ** 2 + (delta[:, 1] / v_scale) ** 2
            weights = 1.0 / np.maximum(distance2, 1.0e-15)
            predicted = float(np.sum(weights * values[mask]) / np.sum(weights))
            errors.append(predicted - values[index])
        return float(np.sqrt(np.mean(np.square(errors))))
