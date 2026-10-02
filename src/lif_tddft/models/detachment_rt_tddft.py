"""RT-TDDFT mean-loss adapter and manuscript sector closure."""

from dataclasses import dataclass

import numpy as np

from lif_tddft.analysis.interpolation import DetachmentProbabilityInterpolator


def detachment_sectors(mean_loss: float) -> tuple[float, float, float]:
    """Map ``Nbar_det`` to (P0, P1, P2) as specified in the detachment-sector model.

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
    sector_interpolators: tuple[DetachmentProbabilityInterpolator, ...]

    @classmethod
    def from_mean_loss_grid(
        cls,
        heights: np.ndarray,
        velocities: np.ndarray,
        mean_loss: np.ndarray,
    ) -> "RTDetachment":
        """Map calculated nodes to probabilities before interpolation.

        Mean electron loss is never interpolated. Only the probabilities of
        zero-, one-, and two-electron loss are interpolated in-domain.
        """

        values = np.asarray(mean_loss, dtype=float)
        if np.any(~np.isfinite(values)) or np.any(values < 0.0) or np.any(values >= 2.0):
            raise ValueError("mean-loss grid must be finite and satisfy 0 <= Nbar_det < 2")
        if values.shape != (len(heights), len(velocities)):
            raise ValueError("mean-loss shape must match the calculated height/velocity nodes")
        sectors = np.asarray([[detachment_sectors(value) for value in row] for row in values])
        return cls(tuple(DetachmentProbabilityInterpolator(heights, velocities, sectors[:, :, i])
                         for i in range(3)))

    def sector_probabilities(self, event: dict) -> tuple[float, float, float]:
        return tuple(interpolator(float(event["surface_height"]), float(event["v_parallel"]))
                     for interpolator in self.sector_interpolators)

    def probability(self, event: dict) -> float:
        p0, _, _ = self.sector_probabilities(event)
        return 1.0 - p0
