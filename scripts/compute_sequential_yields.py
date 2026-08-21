#!/usr/bin/env python3
"""Propagate charge states using a mean-loss grid extracted from RT-TDDFT."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from lif_tddft.models.demkov_parameters import load_demkov_parameters
from lif_tddft.models.charge_state import propagate_ordered_events
from lif_tddft.octopus_results import load_mean_loss_grid, sha256_file


REQUIRED = {
    "trajectory_id", "event_index", "event_time", "surface_height",
    "v_parallel", "energy_defect_au",
}


def read_events(path: Path) -> dict[str, list[dict]]:
    grouped: dict[str, list[dict]] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        missing = REQUIRED - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"encounter table lacks fields: {sorted(missing)}")
        for row in reader:
            event = {
                "trajectory_id": row["trajectory_id"],
                "event_index": int(row["event_index"]),
                "event_time": float(row["event_time"]),
                "surface_height": float(row["surface_height"]),
                "v_parallel": float(row["v_parallel"]),
                "energy_defect_au": float(row["energy_defect_au"]),
            }
            grouped.setdefault(event["trajectory_id"], []).append(event)
    if not grouped:
        raise ValueError("encounter table is empty")
    return grouped


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("encounters", type=Path)
    parser.add_argument("--config", type=Path, default=ROOT / "config" / "demkov_parameters.yaml")
    parser.add_argument("--tddft-grid", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)

    demkov_parameters = load_demkov_parameters(args.config)
    capture = demkov_parameters.capture_model()
    detachment, grid_payload = load_mean_loss_grid(args.tddft_grid)
    summaries = []
    histories = {}
    for trajectory_id, events in read_events(args.encounters).items():
        result = propagate_ordered_events(
            events,
            None,
            capture.probability,
            initial_state="F_neutral",
            detachment_sector_probabilities=detachment.sector_probabilities,
        )
        summaries.append({
            "trajectory_id": trajectory_id,
            "event_count": len(events),
            "final_f_minus_fraction": result["yield_final_hybrid"],
            "final_states": result["states"],
            "transient_f_plus_probability_created": result[
                "transient_f_plus_probability_created"
            ],
        })
        histories[trajectory_id] = result["history"]
    payload = {
        "schema_version": 1,
        "equations": ["manuscript Sec. 2.6", "main Eq. 20", "SI Eq. S14"],
        "data_policy": {
            "normalization": "none",
            "scale_factor": 1.0,
            "value_transform": "none",
            "velocity_specific_adjustments": "none",
        },
        "config": str(args.config),
        "demkov_parameter_definition": {
            "gamma_equation": "(sqrt(2*Et)+sqrt(2*Ep))/2",
            "Et_hartree": demkov_parameters.target_binding_energy_hartree,
            "Ep_hartree": demkov_parameters.projectile_binding_energy_hartree,
            "gamma_bohr_inverse": demkov_parameters.gamma_bohr_inverse,
        },
        "tddft_source": {
            "path": str(args.tddft_grid),
            "sha256": sha256_file(args.tddft_grid),
            "source_kind": grid_payload["source_kind"],
            "equation": grid_payload.get("equation"),
        },
        "summaries": summaries,
        "histories": histories,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
