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
from lif_tddft.models.capture_demkov import DemkovCapture


def gamma_from_config() -> float:
    text = (ROOT / "config" / "demkov_parameters.yaml").read_text(encoding="utf-8")
    line = next(row for row in text.splitlines() if row.strip().startswith("gamma_bohr_inverse:"))
    return float(line.split(":", 1)[1].strip())


def figure4_values() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    wide_path = ROOT / "data" / "processed" / "manuscript" / "detachment_mean_loss_nodes.csv"
    long_path = ROOT / "data" / "processed" / "manuscript" / "figure4_pdet.csv"
    with wide_path.open(encoding="utf-8", newline="") as handle:
        wide = list(csv.DictReader(handle))
    with long_path.open(encoding="utf-8", newline="") as handle:
        long = list(csv.DictReader(handle))
    columns = ["v0p10", "v0p15", "v0p20", "v0p30", "v0p40", "v0p50"]
    source = np.asarray([float(row[column]) for row in wide for column in columns])
    exported = np.asarray([float(row["probability_p_det"]) for row in long])
    source_v010 = np.asarray([float(row["v0p10"]) for row in wide])
    exported_v010 = np.asarray([
        float(row["probability_p_det"])
        for row in long
        if float(row["velocity_au"]) == 0.10
    ])
    return source, exported, source_v010, exported_v010


def capture_v010_error() -> tuple[int, float]:
    path = ROOT / "results" / "tables" / "demkov_s14_event_audit.csv"
    model = DemkovCapture(gamma_from_config())
    differences = []
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            if float(row["velocity_au"]) != 0.10:
                continue
            expected = model.probability({
                "surface_height": float(row["height_bohr"]),
                "v_parallel": 0.10,
                "energy_defect_au": float(row["energy_defect_au"]),
            })
            differences.append(expected - float(row["si_s14_p_capture"]))
    if not differences:
        raise ValueError("the capture audit has no v=0.10 rows")
    return len(differences), float(np.max(np.abs(differences)))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "results" / "tables" / "no_scaling_audit.json",
    )
    args = parser.parse_args(argv)

    policy = json.loads((ROOT / "config" / "data_policy.json").read_text(encoding="utf-8"))[
        "probability_data"
    ]
    assert_no_scaling_policy(policy)
    source, exported, source_v010, exported_v010 = figure4_values()
    validate_probability_values(source, label="Figure 4 Pdet")
    assert_values_identical(source, exported, label="Figure 4 Pdet")
    assert_values_identical(source_v010, exported_v010, label="Figure 4 v=0.10 Pdet")
    capture_count, capture_error = capture_v010_error()
    if capture_error > 1.0e-15:
        raise ValueError("v=0.10 Pcap is not a direct SI Eq. (S14) evaluation")

    payload = {
        "schema_version": 1,
        "passed": True,
        "policy": policy,
        "figure4_probability_values": int(source.size),
        "figure4_v0p10_values": int(source_v010.size),
        "figure4_export_bitwise_identical": True,
        "figure4_v0p10_bitwise_identical": True,
        "capture_v0p10_events_checked": capture_count,
        "capture_v0p10_max_abs_s14_recalculation_error": capture_error,
        "note": "Current manuscript Figure 4 is Pdet; Pcap is Figure 7(a).",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
