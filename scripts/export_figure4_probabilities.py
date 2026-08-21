#!/usr/bin/env python3
"""Export Figure 4 Pdet directly from the production RT-TDDFT mean-loss grid."""

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
from lif_tddft.models.detachment_rt_tddft import detachment_sectors
from lif_tddft.octopus_results import load_mean_loss_grid, sha256_file


POLICY = json.loads((ROOT / "config" / "data_policy.json").read_text(encoding="utf-8"))[
    "probability_data"
]


def pdet_from_mean_loss(value: float) -> float:
    p0, _, _ = detachment_sectors(float(value))
    return 1.0 - p0


def export(source: Path, output: Path, metadata_output: Path) -> None:
    assert_no_scaling_policy(POLICY)
    _, payload = load_mean_loss_grid(source)
    assert_no_scaling_policy(payload["data_policy"])
    heights = np.asarray(payload["heights_bohr"], dtype=float)
    velocities = np.asarray(payload["velocities_au"], dtype=float)
    mean_loss = np.asarray(payload["mean_loss"], dtype=float)

    output_rows = []
    computed_values = []
    exported_values = []
    for height_index, height in enumerate(heights):
        for velocity_index, velocity in enumerate(velocities):
            value = pdet_from_mean_loss(mean_loss[height_index, velocity_index])
            computed_values.append(value)
            output_rows.append({
                "height_bohr": format(float(height), ".17g"),
                "velocity_au": format(float(velocity), ".17g"),
                "probability_p_det": format(value, ".17g"),
                "node_status": "direct_octopus_tddft",
            })
            exported_values.append(float(output_rows[-1]["probability_p_det"]))
    validate_probability_values(np.asarray(computed_values), label="Figure 4 Pdet")
    assert_values_identical(
        np.asarray(computed_values), np.asarray(exported_values), label="Figure 4 Pdet"
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]))
        writer.writeheader()
        writer.writerows(output_rows)
    try:
        source_label = source.resolve().relative_to(ROOT.resolve()).as_posix()
    except ValueError:
        source_label = str(source.resolve())
    metadata = {
        "schema_version": 1,
        "figure": "Figure 4",
        "quantity": "single-collision detachment probability Pdet",
        "source": source_label,
        "source_kind": payload["source_kind"],
        "source_sha256": sha256_file(source),
        "mapping": "Pdet = 1 - P0 from manuscript Sec. 2.6",
        **POLICY,
        "v0p10_special_scaling": False,
    }
    metadata_output.parent.mkdir(parents=True, exist_ok=True)
    metadata_output.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--tddft-grid",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "data" / "processed" / "manuscript" / "figure4_pdet.csv",
    )
    parser.add_argument(
        "--metadata-output",
        type=Path,
        default=ROOT / "data" / "processed" / "manuscript" / "figure4_pdet.metadata.json",
    )
    args = parser.parse_args(argv)
    export(args.tddft_grid, args.output, args.metadata_output)
    print(f"wrote {args.output} without scaling")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
