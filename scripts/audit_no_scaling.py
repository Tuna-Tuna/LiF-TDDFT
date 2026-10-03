#!/usr/bin/env python3
"""Audit identity transforms for manuscript probability data."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from lif_tddft.analysis.data_contract import (
    assert_no_scaling_policy,
    assert_values_identical,
    validate_probability_values,
)
from lif_tddft.models.demkov_parameters import load_demkov_parameters
from lif_tddft.models.detachment_rt_tddft import detachment_probability
from lif_tddft.octopus_results import load_mean_loss_grid


def detachment_values(
    grid_path: Path, long_path: Path
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    _, payload = load_mean_loss_grid(grid_path)
    with long_path.open(encoding="utf-8", newline="") as handle:
        long = list(csv.DictReader(handle))
    velocities = np.asarray(payload["velocities_au"], dtype=float)
    mean_loss = np.asarray(payload["mean_loss"], dtype=float)
    source = np.asarray(
        [detachment_probability(value) for row in mean_loss for value in row],
        dtype=float,
    )
    exported = np.asarray([float(row["probability_p_det"]) for row in long])
    velocity_mask = np.isclose(velocities, 0.10, rtol=0.0, atol=1.0e-12)
    if np.count_nonzero(velocity_mask) != 1:
        raise ValueError("the TDDFT grid must contain exactly one v=0.10 column")
    source_v010 = np.asarray(mean_loss[:, velocity_mask].reshape(-1), dtype=float)
    source_v010 = np.asarray([detachment_probability(value) for value in source_v010])
    exported_v010 = np.asarray([
        float(row["probability_p_det"])
        for row in long
        if float(row["velocity_au"]) == 0.10
    ])
    return source, exported, source_v010, exported_v010


def capture_error(path: Path) -> tuple[int, float]:
    model = load_demkov_parameters(ROOT / "config" / "demkov_parameters.yaml").capture_model()
    differences = []
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            expected = model.probability({
                "surface_height": float(row["height_bohr"]),
                "v_parallel": float(row["v_parallel"]),
                "energy_defect_au": float(row["energy_defect_au"]),
            })
            differences.append(expected - float(row["p_capture"]))
    if not differences:
        raise ValueError("capture table is empty")
    return len(differences), float(np.max(np.abs(differences)))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--tddft-grid",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--detachment-table",
        type=Path,
        default=ROOT / "data" / "processed" / "manuscript" / "detachment_probabilities.csv",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "results" / "tables" / "no_scaling_audit.json",
    )
    parser.add_argument("--capture-table", type=Path)
    args = parser.parse_args(argv)

    policy = json.loads((ROOT / "config" / "data_policy.json").read_text(encoding="utf-8"))[
        "probability_data"
    ]
    assert_no_scaling_policy(policy)
    source, exported, source_v010, exported_v010 = detachment_values(
        args.tddft_grid, args.detachment_table
    )
    validate_probability_values(source, label="Detachment probability Pdet")
    assert_values_identical(source, exported, label="Detachment probability Pdet")
    assert_values_identical(source_v010, exported_v010, label="Detachment v=0.10 Pdet")
    capture_count, capture_difference = 0, None
    if args.capture_table:
        capture_count, capture_difference = capture_error(args.capture_table)
        if not np.isfinite(capture_difference) or capture_difference > 1e-15:
            raise ValueError("Pcap differs from direct evaluation of the capture formula")

    payload = {
        "schema_version": 1,
        "passed": True,
        "policy": policy,
        "detachment_probability_values": int(source.size),
        "detachment_v0p10_values": int(source_v010.size),
        "detachment_export_bitwise_identical": True,
        "detachment_v0p10_bitwise_identical": True,
        "capture_rows_checked": capture_count,
        "capture_max_absolute_error": capture_difference,
        "note": "Detachment nodes use Pdet=Ndet; capture values use direct formula evaluation.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
