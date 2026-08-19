"""RT-TDDFT mean-loss adapter and manuscript sector closure."""

from dataclasses import dataclass

import numpy as np

from lif_tddft.analysis.interpolation import BoundedGridInterpolator


def detachment_sectors(mean_loss: float) -> tuple[float, float, float]:
    """Map ``Nbar_det`` to (P0, P1, P2) as specified in manuscript Sec. 2.6.

    The closure is defined only for ``0 <= Nbar_det < 2``. ``P_det`` for
    loss of the negative-ion state is subsequently ``1 - P0``; it must not
    be confused with the mean number of detached electrons when Nbar_det>1.
    """

    value = float(mean_loss)
    if not 0.0 <= value < 2.0:
        raise ValueError("the manuscript detachment closure requires 0 <= Nbar_det < 2")
    if value <= 1.0:
        return 1.0 - value, value, 0.0
    return 0.0, 2.0 - value, value - 1.0


@dataclass
class RTDetachment:
    interpolator: BoundedGridInterpolator

    @classmethod
    def from_mean_loss_grid(
        cls,
        heights: np.ndarray,
        velocities: np.ndarray,
        mean_loss: np.ndarray,
    ) -> "RTDetachment":
        """Build an unclipped interpolator for ``0 <= Nbar_det < 2`` data."""

        values = np.asarray(mean_loss, dtype=float)
        if np.any(~np.isfinite(values)) or np.any(values < 0.0) or np.any(values >= 2.0):
            raise ValueError("mean-loss grid must be finite and satisfy 0 <= Nbar_det < 2")
        return cls(BoundedGridInterpolator(heights, velocities, values, value_bounds=None))

    def mean_loss(self, event: dict) -> float:
        return self.interpolator(float(event["surface_height"]), float(event["v_parallel"]))

    def sector_probabilities(self, event: dict) -> tuple[float, float, float]:
        return detachment_sectors(self.mean_loss(event))

    def probability(self, event: dict) -> float:
        p0, _, _ = self.sector_probabilities(event)
        return 1.0 - p0
