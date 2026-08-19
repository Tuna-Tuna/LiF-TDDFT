#!/usr/bin/env python3
"""Validate a processed campaign index without inventing missing data."""

import argparse
import csv
import json
from pathlib import Path

REQUIRED = {"run_id", "population_plateau", "flux_closure_relative_rms", "peak_force", "impulse"}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("index", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    with args.index.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        missing = REQUIRED - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"processed index lacks required columns: {sorted(missing)}")
        rows = list(reader)
    if not rows:
        raise ValueError("processed index is empty; missing runs are not filled")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps({"validated_rows": len(rows)}, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
