#!/usr/bin/env python3
"""Compare archived event arrays with the current manuscript/SI equations."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from lif_tddft.models.demkov_parameters import load_demkov_parameters
from lif_tddft.models.charge_state import propagate_ordered_events


def load_csv(root: Path, name: str) -> np.ndarray:
    return np.asarray(np.loadtxt(root / name, delimiter=","), dtype=float)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--legacy-dir", type=Path, default=ROOT / "data" / "processed" / "run_process"
    )
    parser.add_argument(
        "--reference-yields",
        type=Path,
        default=ROOT / "data" / "processed" / "manuscript" / "table_s25_reference_yields.csv",
    )
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results" / "tables")
    args = parser.parse_args(argv)

    demkov_parameters = load_demkov_parameters(
        ROOT / "config" / "demkov_parameters.yaml"
    )
    gamma = demkov_parameters.gamma_bohr_inverse
    demkov = demkov_parameters.capture_model()
    velocity = load_csv(args.legacy_dir, "vcal.csv").reshape(-1)
    heights = load_csv(args.legacy_dir, "Zhr.csv")
    archived_capture = load_csv(args.legacy_dir, "zPcapture.csv")
    energy_defect = load_csv(args.legacy_dir, "zdE.csv")
    archived_detachment = load_csv(args.legacy_dir, "zPloss.csv")
    archived_yields = load_csv(args.legacy_dir, "last_non_zero.csv").reshape(-1)
    reference = np.genfromtxt(args.reference_yields, delimiter=",", names=True, dtype=None, encoding="utf-8")

    rows: list[dict] = []
    summaries = []
    exact_yields = []
    for index, speed in enumerate(velocity):
        columns = np.flatnonzero(np.isfinite(heights[index]) & (heights[index] != 0.0))
        events = []
        errors = []
        for event_index, column in enumerate(columns):
            event = {
                "event_index": event_index,
                "event_time": float(event_index),
                "surface_height": float(heights[index, column]),
                "v_parallel": float(speed),
                "energy_defect_au": float(energy_defect[index, column]),
                "_p_det": float(archived_detachment[index, column]),
            }
            exact = demkov.probability(event)
            archived = float(archived_capture[index, column])
            errors.append(exact - archived)
            event["_p_cap"] = exact
            events.append(event)
            rows.append({
                "velocity_au": speed,
                "event_index": event_index,
                "height_bohr": event["surface_height"],
                "energy_defect_au": event["energy_defect_au"],
                "archived_p_capture": archived,
                "si_s14_p_capture": exact,
                "difference": exact - archived,
            })
        propagated = propagate_ordered_events(
            events,
            lambda event: event["_p_det"],
            lambda event: event["_p_cap"],
            initial_state="F_neutral",
        )
        exact_yield = float(propagated["yield_final_hybrid"])
        exact_yields.append(exact_yield)
        error = np.asarray(errors)
        summaries.append({
            "velocity_au": float(speed),
            "event_count": len(events),
            "max_abs_capture_difference": float(np.max(np.abs(error))),
            "rms_capture_difference": float(np.sqrt(np.mean(error**2))),
            "archived_final_fraction": float(archived_yields[index]),
            "si_s14_with_archived_detachment_final_fraction": exact_yield,
            "table_s25_reference_fraction": float(reference["final_f_minus_fraction"][index]),
        })

    args.output_dir.mkdir(parents=True, exist_ok=True)
    event_path = args.output_dir / "demkov_s14_event_audit.csv"
    with event_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    summary_path = args.output_dir / "manuscript_formula_audit.csv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summaries[0]))
        writer.writeheader()
        writer.writerows(summaries)
    maximum = max(row["max_abs_capture_difference"] for row in summaries)
    payload = {
        "schema_version": 1,
        "production_formula": "SI Eq. S14 without an additional height envelope",
        "gamma_formula": "(sqrt(2*Et)+sqrt(2*Ep))/2",
        "Et_hartree": demkov_parameters.target_binding_energy_hartree,
        "Ep_hartree": demkov_parameters.projectile_binding_energy_hartree,
        "gamma_bohr_inverse": gamma,
        "maximum_archived_vs_s14_event_probability_difference": maximum,
        "archived_capture_array_matches_si_s14": maximum <= 1.0e-12,
        "public_release_blocker": maximum > 1.0e-12,
        "reason": (
            "The archived zPcapture array contains an additional height-dependent factor that "
            "is absent from SI Eq. S14. Table S25 reference yields also use a later detachment "
            "surface not available as encounter-resolved release data."
        ),
        "outputs": [
            event_path.relative_to(ROOT).as_posix(),
            summary_path.relative_to(ROOT).as_posix(),
        ],
    }
    (args.output_dir / "manuscript_formula_audit.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
