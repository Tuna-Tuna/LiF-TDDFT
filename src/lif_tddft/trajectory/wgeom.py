"""Pure geometric access weights; never a final-yield integral."""

from __future__ import annotations

import numpy as np


def wgeom(events: list[dict], height_edges: np.ndarray, velocity_edges: np.ndarray) -> dict[str, object]:
    if not events:
        raise ValueError("at least one event is required")
    heights = np.asarray([row["surface_height"] for row in events], dtype=float)
    velocities = np.asarray([row["v_parallel"] for row in events], dtype=float)
    weights = np.asarray([row.get("trajectory_weight", 1.0) for row in events], dtype=float)
    trajectory_weights: dict[str, float] = {}
    for row in events:
        trajectory_weights[str(row["trajectory_id"])] = float(row.get("trajectory_weight", 1.0))
    incident_weight = sum(trajectory_weights.values())
    if incident_weight <= 0:
        raise ValueError("incident trajectory weight must be positive")
    counts, _, _ = np.histogram2d(heights, velocities, bins=(height_edges, velocity_edges), weights=weights)
    density = counts / incident_weight
    thresholds = {str(limit): float(weights[heights < limit].sum() / incident_weight) for limit in (3.5, 3.0, 2.5, 2.0)}
    ratio = np.asarray([abs(row["v_perpendicular"]) / max(abs(row["v_parallel"]), 1.0e-15) for row in events])
    return {
        "event_density": density,
        "mean_effective_collision_count": float(weights.sum() / incident_weight),
        "access_fraction_by_height": thresholds,
        "v_perpendicular_over_v_parallel": ratio,
    }


def database_coverage(events: list[dict], h_range: tuple[float, float], v_range: tuple[float, float]) -> dict[str, float]:
    weights = np.asarray([row.get("trajectory_weight", 1.0) for row in events], dtype=float)
    inside = np.asarray([
        h_range[0] <= row["surface_height"] <= h_range[1]
        and v_range[0] <= row["v_parallel"] <= v_range[1]
        for row in events
    ])
    total = float(weights.sum())
    covered = float(weights[inside].sum())
    return {"covered_weight": covered, "out_of_domain_weight": total - covered, "coverage_fraction": covered / total if total else 0.0}
