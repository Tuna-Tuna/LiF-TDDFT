#!/usr/bin/env python3
"""Build the production mean-loss grid directly from paired Octopus outputs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from lif_tddft.analysis.data_contract import assert_no_scaling_policy
from lif_tddft.octopus_results import extract_paired_plateau, load_completed_run, sha256_file


DATA_POLICY = {
    "normalization": "none",
    "scale_factor": 1.0,
    "value_transform": "none",
    "velocity_specific_adjustments": "none",
}


def pair_runs(campaign_root: Path) -> list[tuple[dict, Path, Path, Path, Path]]:
    grouped: dict[tuple[float, float, str], dict[str, tuple[dict, Path, Path]]] = {}
    for manifest_path in sorted(campaign_root.glob("*/manifest.json")):
        run_dir = manifest_path.parent
        manifest, stage = load_completed_run(run_dir)
        run = manifest["run"]
        key = (
            float(run["height_bohr"]),
            float(run["relative_velocity_au"]),
            str(run["cluster_model"]),
        )
        variant = str(run["variant"])
        if variant in grouped.setdefault(key, {}):
            raise ValueError(f"duplicate {variant} run for {key}")
        grouped[key][variant] = (manifest, stage, manifest_path)
    if not grouped:
        raise ValueError(f"no completed run manifests found below {campaign_root}")
    pairs = []
    for key, variants in sorted(grouped.items()):
        if set(variants) != {"production", "isolated_projectile"}:
            raise ValueError(f"paired TDDFT runs are incomplete for {key}: {sorted(variants)}")
        production, isolated = variants["production"], variants["isolated_projectile"]
        pairs.append((production[0], production[1], isolated[1], production[2], isolated[2]))
    return pairs


def build_grid(
    campaign_root: Path,
    *,
    radius_angstrom: float,
    tail_fraction: float,
    transition_half_width_angstrom: float = 0.3,
) -> dict[str, object]:
    if radius_angstrom <= 0.0:
        raise ValueError("radius_angstrom must be positive")
    if not 0.0 < tail_fraction <= 1.0:
        raise ValueError("tail_fraction must lie in (0,1]")
    assert_no_scaling_policy(DATA_POLICY)
    extracted = []
    for manifest, production_stage, isolated_stage, production_manifest, isolated_manifest in pair_runs(campaign_root):
        run = manifest["run"]
        mean_loss, provenance = extract_paired_plateau(
            production_stage,
            isolated_stage,
            radius=radius_angstrom,
            transition_half_width=transition_half_width_angstrom,
            tail_fraction=tail_fraction,
        )
        extracted.append(
            {
                "height_bohr": float(run["height_bohr"]),
                "velocity_au": float(run["relative_velocity_au"]),
                "cluster_model": str(run["cluster_model"]),
                "mean_loss": mean_loss,
                "production_manifest_sha256": sha256_file(production_manifest),
                "isolated_manifest_sha256": sha256_file(isolated_manifest),
                "provenance": provenance,
            }
        )
    models = sorted({item["cluster_model"] for item in extracted})
    if len(models) != 1:
        raise ValueError(f"one extraction output must contain one cluster model; found {models}")
    heights = sorted({float(item["height_bohr"]) for item in extracted})
    velocities = sorted({float(item["velocity_au"]) for item in extracted})
    lookup = {(float(item["height_bohr"]), float(item["velocity_au"])): item for item in extracted}
    expected = {(height, velocity) for height in heights for velocity in velocities}
    if set(lookup) != expected:
        raise ValueError("completed paired runs do not form a rectangular height/velocity grid")
    values = np.asarray(
        [[float(lookup[(height, velocity)]["mean_loss"]) for velocity in velocities] for height in heights]
    )
    if np.any(~np.isfinite(values)) or np.any(values < 0.0) or np.any(values > 1.0):
        raise ValueError("extracted mean loss must be finite and satisfy 0 <= Ndet <= 1")
    return {
        "schema_version": 1,
        "source_kind": "octopus_rt_tddft_density",
        "equation": "paired interacting-minus-isolated population change",
        "cluster_model": models[0],
        "bound_radius_angstrom": radius_angstrom,
        "transition_half_width_angstrom": transition_half_width_angstrom,
        "tail_fraction": tail_fraction,
        "data_policy": DATA_POLICY,
        "heights_bohr": heights,
        "velocities_au": velocities,
        "mean_loss": values.tolist(),
        "runs": extracted,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("campaign_root", type=Path)
    parser.add_argument("--radius-angstrom", type=float, default=3.5)
    parser.add_argument("--transition-half-width-angstrom", type=float, default=0.3)
    parser.add_argument("--tail-fraction", type=float, default=0.2)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    payload = build_grid(
        args.campaign_root,
        radius_angstrom=args.radius_angstrom,
        transition_half_width_angstrom=args.transition_half_width_angstrom,
        tail_fraction=args.tail_fraction,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
