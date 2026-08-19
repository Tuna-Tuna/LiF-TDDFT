#!/usr/bin/env python3
"""Regenerate the tracked one-input-per-height Octopus examples."""

from __future__ import annotations

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from lif_tddft.campaign import expand_campaign, load_config
from lif_tddft.octopus_input import render_input


def token(value: float) -> str:
    return f"{value:.1f}".replace(".", "_")


def main() -> int:
    config = load_config(ROOT / "config" / "campaign_revision.yaml")
    specs = [
        item for item in expand_campaign(config)
        if item.variant == "production" and item.relative_velocity_au == 0.30
    ]
    written = []
    for spec in specs:
        target = ROOT / "octopus" / "inputs" / f"h{token(spec.height_bohr)}" / "inp"
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(
            render_input(
                ROOT / "octopus" / "templates" / "single_flyby.inp.j2",
                spec,
                config,
                move_ions=True,
                constant_velocity=False,
                calculation_mode="td",
            ),
            encoding="utf-8",
        )
        written.append(target.relative_to(ROOT).as_posix())
    (ROOT / "octopus" / "inputs" / "manifest.json").write_text(
        json.dumps({
            "schema_version": 1,
            "purpose": "tracked v=0.30 examples; full 96-run matrix is generated from config",
            "files": written,
        }, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"wrote {len(written)} reference inputs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
