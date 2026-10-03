"""Single-electron detachment: Pdet equals the corrected outgoing loss."""
from dataclasses import dataclass
import numpy as np
from lif_tddft.analysis.interpolation import DetachmentProbabilityInterpolator


def detachment_probability(mean_loss: float) -> float:
    """PRA Sec. II D: only zero or one electron can be removed."""
    value = float(mean_loss)
    if not np.isfinite(value) or not 0.0 <= value <= 1.0:
        raise ValueError("single-electron detachment requires finite 0 <= Ndet <= 1; clipping is forbidden")
    return value


@dataclass
class RTDetachment:
    interpolator: DetachmentProbabilityInterpolator

    @classmethod
    def from_mean_loss_grid(cls, heights, velocities, mean_loss):
        values = np.asarray(mean_loss, dtype=float)
        if np.any(~np.isfinite(values)) or np.any(values < 0.0) or np.any(values > 1.0):
            raise ValueError("single-electron grid must be finite and satisfy 0 <= Ndet <= 1")
        return cls(DetachmentProbabilityInterpolator(heights, velocities, values))

    def probability(self, event: dict) -> float:
        return self.interpolator(float(event['surface_height']), float(event['v_parallel']))
