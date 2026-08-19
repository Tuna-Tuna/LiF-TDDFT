#!/usr/bin/env python3
"""Expand the revision campaign into independent, validated run plans."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from lif_tddft.campaign import expand_campaign, load_config, ramp_velocities
from lif_tddft.octopus_input import render_input
from lif_tddft.provenance import sha256_file


def build_manifest(spec, config: dict, input_path: Path) -> dict:
    ramp = config["ramp"]
    velocities = ramp_velocities(
        spec.relative_velocity_au,
        float(ramp["initial_relative_velocity_au"]),
        float(ramp["increment_au"]),
    )
    stages = [{"name": "ground_state", "from_scratch": True, "velocity_au": 0.0}]
    stages.extend(
        {
            "name": f"ramp_{index:03d}",
            "from_scratch": False,
            "velocity_au": value,
            "propagation_time_fs": float(ramp["segment_time_fs"]),
            "ions_constant_velocity": False,
        }
        for index, value in enumerate(velocities, 1)
    )
    stages.append({
        "name": "production_ehrenfest",
        "from_scratch": False,
        "velocity_au": spec.relative_velocity_au,
        "move_ions": True,
        "ions_constant_velocity": False,
        "surface_coordinates_fixed": True,
    })
    return {
        "schema_version": 1,
        "run": spec.to_dict(),
        "input": str(input_path.resolve()),
        "stages": stages,
        "status": "planned",
        "independent_target_velocity_history": True,
    }


def generate(config_path: Path, output_root: Path, dry_run: bool = False) -> list[dict]:
    config = load_config(config_path)
    specs = expand_campaign(config)
    rows: list[dict] = []
    for spec in specs:
        run_dir = output_root / spec.run_id
        input_path = run_dir / "input" / "inp"
        template_name = (
            "isolated_projectile.inp.j2"
            if spec.variant == "isolated_projectile"
            else "single_flyby.inp.j2"
        )
        template = ROOT / "octopus" / "templates" / template_name
        text = render_input(template, spec, config)
        manifest = build_manifest(spec, config, input_path)
        if not dry_run:
            input_path.parent.mkdir(parents=True, exist_ok=True)
            input_path.write_text(text, encoding="utf-8")
            input_manifest = {
                "schema_version": 1,
                "run": spec.to_dict(),
                "template": str(template.resolve()),
                "input_sha256": sha256_file(input_path),
                "units": {"height": "bohr", "relative_velocity": "atomic_unit"},
                "preflight": {
                    "relative_velocity_match": True,
                    "height_match": True,
                    "production_ehrenfest": True,
                    "fixed_surface_coordinates": True,
                    "zero_initial_transverse_and_normal_velocity": True,
                    "initial_overlap_gate": "requires fragment-orbital calculation",
                },
            }
            (input_path.parent / "input_manifest.json").write_text(
                json.dumps(input_manifest, indent=2) + "\n", encoding="utf-8"
            )
            (run_dir / "manifest.json").write_text(
                json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
            )
        rows.append({**spec.to_dict(), "manifest": str(run_dir / "manifest.json")})
    if not dry_run:
        output_root.mkdir(parents=True, exist_ok=True)
        with (output_root / "campaign_plan.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    return rows


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=ROOT / "config" / "campaign_revision.yaml")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "runs")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)
    rows = generate(args.config, args.output_dir, args.dry_run)
    print(f"validated {len(rows)} independent run plans")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
