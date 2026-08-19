#!/usr/bin/env python3
"""Revision figure gate: figures consume frozen tables only."""

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from lif_tddft.analysis.data_contract import assert_no_scaling_policy


def load_metadata(path: Path) -> dict:
    text = path.read_text(encoding="utf-8")
    try:
        value = json.loads(text)
    except json.JSONDecodeError:
        try:
            import yaml
        except ImportError as exc:
            raise RuntimeError("YAML figure metadata requires PyYAML") from exc
        value = yaml.safe_load(text)
    if not isinstance(value, dict):
        raise ValueError("figure metadata must be a mapping")
    return value


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("table", type=Path)
    parser.add_argument("--metadata", type=Path, required=True)
    args = parser.parse_args(argv)
    if not args.table.is_file() or args.table.stat().st_size == 0:
        raise FileNotFoundError("a non-empty machine-readable source table is required")
    if not args.metadata.is_file() or args.metadata.stat().st_size == 0:
        raise FileNotFoundError("figure metadata is required")
    metadata = load_metadata(args.metadata)
    policy = metadata.get("probability_data", metadata)
    if not isinstance(policy, dict):
        raise ValueError("probability_data policy must be a mapping")
    assert_no_scaling_policy(policy)
    print("inputs validated with identity probability transform; plotting is project-specific")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
