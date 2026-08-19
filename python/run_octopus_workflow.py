#!/usr/bin/env python3
"""Fail-fast ground-state, ramp, and Ehrenfest production runner."""

from __future__ import annotations

import argparse
import json
import os
import shlex
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from lif_tddft.campaign import RunSpec, load_config
from lif_tddft.octopus_input import render_input
from lif_tddft.provenance import build_provenance, software_version, write_json
from lif_tddft.restart import validate_restart


def _status(path: Path, state: str, stages: list[dict], error: str | None = None) -> None:
    write_json(
        path,
        {
            "state": state,
            "updated_utc": datetime.now(timezone.utc).isoformat(),
            "stages": stages,
            "error": error,
        },
    )


def run_plan(
    plan_path: Path,
    config_path: Path,
    octopus: str,
    mpi_command: list[str],
    *,
    dry_run: bool = False,
) -> None:
    manifest = json.loads(plan_path.read_text(encoding="utf-8"))
    config = load_config(config_path)
    spec = RunSpec(**manifest["run"])
    run_dir = plan_path.parent
    template_name = (
        "isolated_projectile.inp.j2"
        if spec.variant == "isolated_projectile"
        else "single_flyby.inp.j2"
    )
    template = ROOT / "octopus" / "templates" / template_name
    status_path = run_dir / "run_status.json"
    completed: list[dict] = []
    _status(status_path, "running" if not dry_run else "validated", completed)
    previous_stage: Path | None = None
    for stage in manifest["stages"]:
        stage_dir = run_dir / "stages" / stage["name"]
        velocity = float(stage["velocity_au"])
        is_ground_state = stage["name"] == "ground_state"
        is_ramp = stage["name"].startswith("ramp_")
        constant = bool(stage.get("ions_constant_velocity", False))
        from_scratch = "yes" if bool(stage.get("from_scratch")) else "no"
        text = render_input(
            template,
            spec,
            config,
            relative_velocity_au=velocity,
            from_scratch=from_scratch,
            move_ions=not is_ground_state,
            constant_velocity=constant,
            calculation_mode="gs" if is_ground_state else "td",
            propagation_time_expression=(
                f"{float(config['ramp']['segment_time_fs'])}*fs" if is_ramp else None
            ),
        )
        if dry_run:
            completed.append({"name": stage["name"], "state": "validated"})
            continue
        stage_dir.mkdir(parents=True, exist_ok=False)
        (stage_dir / "inp").write_text(text, encoding="utf-8")
        required_pseudos = ["F.oncvpsp.psp8"]
        if spec.variant != "isolated_projectile":
            required_pseudos.append("Li.oncvpsp.psp8")
        for pseudo_name in required_pseudos:
            pseudo_source = run_dir / "input" / pseudo_name
            if not pseudo_source.is_file():
                _status(status_path, "failed", completed, f"missing pseudopotential: {pseudo_source}")
                raise FileNotFoundError(pseudo_source)
            shutil.copy2(pseudo_source, stage_dir / pseudo_name)
        if previous_stage is not None:
            source = previous_stage / "restart"
            if not source.is_dir():
                raise RuntimeError(f"required restart is absent: {source}")
            shutil.copytree(source, stage_dir / "restart")
        command = [*mpi_command, octopus]
        with (stage_dir / "run.log").open("w", encoding="utf-8") as log:
            result = subprocess.run(command, cwd=stage_dir, stdout=log, stderr=subprocess.STDOUT)
        if result.returncode != 0:
            _status(status_path, "failed", completed, f"{stage['name']} exit={result.returncode}")
            raise RuntimeError(f"Octopus failed in {stage['name']}")
        check = validate_restart(
            stage_dir,
            require_coordinates=stage["name"] != "ground_state",
        )
        if not check["valid"]:
            _status(status_path, "failed", completed, "; ".join(check["issues"]))
            raise RuntimeError(f"incomplete restart after {stage['name']}: {check['issues']}")
        completed.append({"name": stage["name"], "state": "complete", "restart": check})
        _status(status_path, "running", completed)
        previous_stage = stage_dir
    if dry_run:
        _status(status_path, "validated", completed)
        return
    pseudo_names = ["F.oncvpsp.psp8"]
    if spec.variant != "isolated_projectile":
        pseudo_names.append("Li.oncvpsp.psp8")
    pseudo = [run_dir / "input" / name for name in pseudo_names]
    provenance = build_provenance(
        ROOT,
        [run_dir / "input" / "inp"],
        pseudo,
        octopus_version=software_version(octopus, "--version"),
        mpi_command=mpi_command,
    )
    write_json(run_dir / "provenance.json", provenance)
    _status(status_path, "complete", completed)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--plan", required=True, type=Path)
    parser.add_argument("--config", type=Path, default=ROOT / "config" / "campaign_revision.yaml")
    parser.add_argument("--octopus", default=os.environ.get("OCTOPUS_PATH", "octopus"))
    parser.add_argument("--mpi-cmd", default=os.environ.get("OCTOPUS_MPI_CMD", ""))
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)
    run_plan(args.plan, args.config, args.octopus, shlex.split(args.mpi_cmd), dry_run=args.dry_run)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
