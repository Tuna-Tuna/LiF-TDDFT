"""Named uncertainty sources; no anonymous aggregate band."""

from __future__ import annotations

import math

ALLOWED_SOURCES = {
    "uncertainty_tddft_numerical",
    "uncertainty_cluster_size",
    "uncertainty_force_interpolation",
    "uncertainty_trajectory_sampling",
    "uncertainty_demkov_model",
}


def combine(components: dict[str, float]) -> dict[str, object]:
    unknown = set(components) - ALLOWED_SOURCES
    if unknown:
        raise ValueError(f"unrecognized uncertainty sources: {sorted(unknown)}")
    if any(value < 0 for value in components.values()):
        raise ValueError("uncertainties must be non-negative")
    return {"components": dict(components), "quadrature_total": math.sqrt(sum(value * value for value in components.values()))}
