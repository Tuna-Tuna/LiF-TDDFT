#!/usr/bin/env python3
"""Evaluate SI Eq. (S14) for a Figure 7(a)-style grid without scaling."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from lif_tddft.analysis.data_contract import validate_probability_values
from lif_tddft.models.capture_demkov import DemkovCapture


REQUIRED = {"height_bohr", "v_parallel", "energy_defect_au"}


def gamma_from_config(path: Path) -> float:
    line = next(
        row for row in path.read_text(encoding="utf-8").splitlines()
        if row.strip().startswith("gamma_bohr_inverse:")
    )
    return float(line.split(":", 1)[1].strip())


def compute(source: Path, output: Path, config: Path) -> None:
    model = DemkovCapture(gamma_from_config(config))
    with source.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        missing = REQUIRED - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"capture grid lacks columns: {sorted(missing)}")
        fieldnames = list(reader.fieldnames or [])
        if "p_capture" in fieldnames:
            raise ValueError("capture input must not contain a precomputed p_capture column")
        rows = list(reader)
    if not rows:
        raise ValueError("capture grid is empty")
    probabilities = []
    for row in rows:
        event = {
            "surface_height": float(row["height_bohr"]),
            "v_parallel": float(row["v_parallel"]),
            "energy_defect_au": float(row["energy_defect_au"]),
        }
        probability = model.probability(event)
        probabilities.append(probability)
        row["p_capture"] = repr(probability)
    validate_probability_values(np.asarray(probabilities), label="Figure 7(a) Pcap")
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=[*fieldnames, "p_capture"])
        writer.writeheader()
        writer.writerows(rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--config", type=Path, default=ROOT / "config" / "demkov_parameters.yaml"
    )
    args = parser.parse_args(argv)
    compute(args.source, args.output, args.config)
    print(f"wrote {args.output} from SI Eq. (S14) without scaling")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
