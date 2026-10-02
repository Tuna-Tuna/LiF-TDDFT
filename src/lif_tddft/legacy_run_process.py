"""Read and audit variables extracted from the legacy run_process.m result."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from .models.charge_state import propagate_ordered_events
from .provenance import sha256_file
from .trajectory.wgeom import database_coverage


REQUIRED_ARRAYS = {
    "vcal", "zcal", "result_hsemt", "result_esemt", "result_loss",
    "Zhr", "zPcapture", "zPloss", "Pfin", "last_non_zero", "Zturn", "expe",
}


def _read_csv(path: Path) -> np.ndarray:
    value = np.loadtxt(path, delimiter=",")
    return np.asarray(value, dtype=float)


def load_extracted(directory: str | Path) -> dict[str, Any]:
    root = Path(directory)
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    arrays = {record["name"]: _read_csv(root / record["file"]) for record in manifest["arrays"]}
    missing = REQUIRED_ARRAYS - arrays.keys()
    if missing:
        raise ValueError(f"extraction is missing required arrays: {sorted(missing)}")
    arrays["vcal"] = arrays["vcal"].reshape(-1)
    arrays["zcal"] = arrays["zcal"].reshape(-1)
    arrays["last_non_zero"] = arrays["last_non_zero"].reshape(-1)
    arrays["Zturn"] = arrays["Zturn"].reshape(-1)
    if arrays["result_loss"].shape != (arrays["zcal"].size, arrays["vcal"].size):
        raise ValueError("result_loss shape is inconsistent with zcal/vcal")
    event_shape = (arrays["vcal"].size, arrays["Zhr"].shape[1])
    for name in ("Zhr", "zPcapture", "zPloss", "Pfin"):
        if arrays[name].shape != event_shape:
            raise ValueError(f"{name} has unexpected shape {arrays[name].shape}")
    return {"manifest": manifest, "arrays": arrays}


def _legacy_events(arrays: dict[str, np.ndarray], velocity_index: int) -> list[dict[str, Any]]:
    heights = arrays["Zhr"][velocity_index]
    valid_indices = np.flatnonzero(np.isfinite(heights) & (heights != 0.0))
    velocity = float(arrays["vcal"][velocity_index])
    events = []
    for event_index, column in enumerate(valid_indices):
        events.append(
            {
                "trajectory_id": f"legacy_v{velocity:g}",
                "trajectory_weight": 1.0,
                "event_index": event_index,
                "event_time": float(event_index),
                "surface_height": float(heights[column]),
                "v_parallel": velocity,
                "v_perpendicular": float("nan"),
                "_p_det": float(arrays["zPloss"][velocity_index, column]),
                "_p_cap": float(arrays["zPcapture"][velocity_index, column]),
            }
        )
    return events


def replay(directory: str | Path) -> dict[str, Any]:
    loaded = load_extracted(directory)
    arrays = loaded["arrays"]
    summaries = []
    all_events: list[dict[str, Any]] = []
    for velocity_index, velocity in enumerate(arrays["vcal"]):
        events = _legacy_events(arrays, velocity_index)
        all_events.extend(events)
        propagated = propagate_ordered_events(
            events,
            lambda event: event["_p_det"],
            lambda event: event["_p_cap"],
            initial_state="F_neutral",
        )
        below_domain = sum(event["surface_height"] < arrays["zcal"][0] for event in events)
        above_domain = sum(event["surface_height"] > arrays["zcal"][-1] for event in events)
        out_of_domain = below_domain + above_domain
        saved_final = float(arrays["last_non_zero"][velocity_index])
        python_final = float(propagated["yield_final_hybrid"])
        summaries.append(
            {
                "velocity_au": float(velocity),
                "event_count": len(events),
                "minimum_height_bohr": min(event["surface_height"] for event in events),
                "maximum_height_bohr": max(event["surface_height"] for event in events),
                "out_of_domain_events": out_of_domain,
                "below_domain_events": below_domain,
                "above_domain_events": above_domain,
                "td_database_coverage_fraction": (len(events) - out_of_domain) / len(events),
                "matlab_final_yield": saved_final,
                "python_replayed_final_yield": python_final,
                "replay_absolute_error": abs(python_final - saved_final),
                "turning_height_bohr": float(arrays["Zturn"][velocity_index]),
            }
        )
    coverage = database_coverage(
        all_events,
        (float(arrays["zcal"][0]), float(arrays["zcal"][-1])),
        (float(arrays["vcal"][0]), float(arrays["vcal"][-1])),
    )
    experiment = np.asarray(arrays["expe"], dtype=float)
    differences = []
    matched_points = 0
    for velocity, observed, *_ in experiment:
        matches = np.flatnonzero(np.isclose(arrays["vcal"], velocity, rtol=0, atol=1e-12))
        if matches.size > 1:
            raise ValueError("ambiguous calculated velocity in experimental comparison")
        if matches.size == 1:
            differences.append(float(arrays["last_non_zero"][matches[0]]) - float(observed))
            matched_points += 1
    experiment_rmse = float(np.sqrt(np.mean(np.square(differences)))) if differences else None
    source_path = Path(loaded["manifest"]["source_mat"])
    return {
        "schema_version": 1,
        "source_mat": str(source_path),
        "source_mat_sha256": sha256_file(source_path) if source_path.is_file() else None,
        "velocity_summaries": summaries,
        "event_weight_coverage": coverage,
        "experiment_points_at_calculated_velocities": matched_points,
        "experiment_points_unmatched": int(len(experiment) - matched_points),
        "unscaled_theory_experiment_rmse": experiment_rmse,
        "limitations": [
            "Saved capture arrays are used only to replay the supplied results.",
            "Legacy result_loss values are audited unchanged; Python does not clip or rescale them.",
            "Saved event arrays do not contain normal velocity; v_perpendicular is unavailable.",
            "Only supplied event probabilities are replayed; no missing values are generated.",
            "Experimental comparison uses matching calculated velocities only; unmatched points are reported.",
        ],
    }
