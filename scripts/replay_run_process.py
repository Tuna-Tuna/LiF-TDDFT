#!/usr/bin/env python3
"""Replay saved run_process probabilities through the revision Python model."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from lif_tddft.legacy_run_process import replay


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("extracted", type=Path)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results" / "tables")
    args = parser.parse_args(argv)
    result = replay(args.extracted)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    json_path = args.output_dir / "run_process_replay.json"
    csv_path = args.output_dir / "run_process_replay_summary.csv"
    json_path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    rows = result["velocity_summaries"]
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({
        "json": str(json_path),
        "csv": str(csv_path),
        "velocities": len(rows),
        "maximum_replay_error": max(row["replay_absolute_error"] for row in rows),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
