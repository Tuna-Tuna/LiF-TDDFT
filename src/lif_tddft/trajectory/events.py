"""Ordered closest-approach event extraction."""

from __future__ import annotations

import numpy as np


EVENT_FIELDS = (
    "trajectory_id", "trajectory_weight", "incident_velocity", "incident_angle",
    "azimuth", "initial_impact_parameter", "event_index", "site_id", "event_time",
    "x", "y", "z", "surface_height", "v_parallel", "v_perpendicular",
    "local_incident_angle", "minimum_lateral_distance",
)


def identify_events(
    trajectory_id: str,
    times: np.ndarray,
    positions: np.ndarray,
    velocities: np.ndarray,
    sites: dict[str, np.ndarray],
    *,
    trajectory_weight: float = 1.0,
    incident_angle: float = 0.0,
    azimuth: float = 0.0,
    initial_impact_parameter: float = 0.0,
    surface_z: float = 0.0,
) -> list[dict[str, float | str | int]]:
    t = np.asarray(times, dtype=float)
    r = np.asarray(positions, dtype=float)
    v = np.asarray(velocities, dtype=float)
    if r.shape != v.shape or r.shape != (t.size, 3):
        raise ValueError("times, positions, and velocities are inconsistent")
    rows = []
    for site_id, coordinate in sites.items():
        site = np.asarray(coordinate, dtype=float).reshape(3)
        lateral2 = np.sum((r[:, :2] - site[:2]) ** 2, axis=1)
        candidates = [i for i in range(1, t.size - 1) if lateral2[i] <= lateral2[i-1] and lateral2[i] < lateral2[i+1]]
        for index in candidates:
            parallel = float(np.linalg.norm(v[index, :2]))
            perpendicular = float(v[index, 2])
            local_angle = float(np.arctan2(abs(perpendicular), max(parallel, 1.0e-15)))
            rows.append({
                "trajectory_id": trajectory_id,
                "trajectory_weight": float(trajectory_weight),
                "incident_velocity": float(np.linalg.norm(v[0])),
                "incident_angle": float(incident_angle),
                "azimuth": float(azimuth),
                "initial_impact_parameter": float(initial_impact_parameter),
                "site_id": site_id,
                "event_time": float(t[index]),
                "x": float(r[index, 0]), "y": float(r[index, 1]), "z": float(r[index, 2]),
                "surface_height": float(r[index, 2] - surface_z),
                "v_parallel": parallel,
                "v_perpendicular": perpendicular,
                "local_incident_angle": local_angle,
                "minimum_lateral_distance": float(np.sqrt(lateral2[index])),
            })
    rows.sort(key=lambda row: (float(row["event_time"]), str(row["site_id"])))
    for event_index, row in enumerate(rows):
        row["event_index"] = event_index
    return rows
