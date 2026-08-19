#!/usr/bin/env python3
"""Export Figure 4 Pdet nodes without changing any probability value."""

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


VELOCITY_COLUMNS = {
    "v0p10": 0.10,
    "v0p15": 0.15,
    "v0p20": 0.20,
    "v0p30": 0.30,
    "v0p40": 0.40,
    "v0p50": 0.50,
}
POLICY = json.loads((ROOT / "config" / "data_policy.json").read_text(encoding="utf-8"))[
    "probability_data"
]


def export(source: Path, output: Path, metadata_output: Path) -> None:
    assert_no_scaling_policy(POLICY)
    with source.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        missing = {"height_bohr", "node_status", *VELOCITY_COLUMNS} - set(
            reader.fieldnames or []
        )
        if missing:
            raise ValueError(f"Figure 4 source lacks columns: {sorted(missing)}")
        source_rows = list(reader)
    if not source_rows:
        raise ValueError("Figure 4 source table is empty")

    output_rows = []
    source_values = []
    exported_values = []
    for row in source_rows:
        for column, velocity in VELOCITY_COLUMNS.items():
            value_text = row[column]
            value = float(value_text)
            source_values.append(value)
            output_rows.append({
                "height_bohr": row["height_bohr"],
                "velocity_au": f"{velocity:.2f}",
                "probability_p_det": value_text,
                "node_status": row["node_status"],
            })
            exported_values.append(float(value_text))
    validate_probability_values(np.asarray(source_values), label="Figure 4 Pdet")
    assert_values_identical(
        np.asarray(source_values), np.asarray(exported_values), label="Figure 4 Pdet"
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
        **POLICY,
        "v0p10_special_scaling": False,
    }
    metadata_output.parent.mkdir(parents=True, exist_ok=True)
    metadata_output.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--source",
        type=Path,
        default=ROOT / "data" / "processed" / "manuscript" / "detachment_mean_loss_nodes.csv",
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
    export(args.source, args.output, args.metadata_output)
    print(f"wrote {args.output} without scaling")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
