"""Strict extraction of projectile populations from Octopus RT-TDDFT output."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import re
from typing import Callable

import numpy as np

from lif_tddft.analysis.population import (
    moving_sphere_population,
    paired_projectile_mean_loss,
    plateau_statistics,
)
from lif_tddft.analysis.data_contract import assert_no_scaling_policy
from lif_tddft.models.detachment_rt_tddft import RTDetachment, detachment_sectors


_TD_DIRECTORY = re.compile(r"^td\.(\d+)$")


@dataclass(frozen=True)
class CoordinateSeries:
    iterations: np.ndarray
    times: np.ndarray
    projectile_positions: np.ndarray


@dataclass(frozen=True)
class PopulationSeries:
    iterations: np.ndarray
    times: np.ndarray
    projectile_positions: np.ndarray
    populations: np.ndarray
    source_files: tuple[Path, ...]


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_td_coordinates(
    path: str | Path,
    *,
    natoms: int,
    projectile_index: int = 0,
) -> CoordinateSeries:
    """Read iteration, time, and projectile position from ``td.general/coordinates``.

    Octopus writes positions first, followed by velocities and forces. Only the
    position block is consumed here. Values remain in ``UnitsOutput`` length
    units, matching the spatial density grid written by the same run.
    """

    if natoms <= 0 or not 0 <= projectile_index < natoms:
        raise ValueError("invalid atom count or projectile index")
    records: list[tuple[int, float, np.ndarray]] = []
    minimum_columns = 2 + 3 * natoms
    for line_number, raw in enumerate(Path(path).read_text(encoding="utf-8").splitlines(), 1):
        stripped = raw.strip()
        if not stripped or stripped.startswith("#"):
            continue
        fields = stripped.split()
        if len(fields) < minimum_columns:
            raise ValueError(
                f"{path}:{line_number} has {len(fields)} columns; expected at least {minimum_columns}"
            )
        try:
            iteration = int(fields[0])
            time = float(fields[1].replace("D", "E").replace("d", "e"))
            positions = np.asarray(
                [float(value.replace("D", "E").replace("d", "e")) for value in fields[2:minimum_columns]],
                dtype=float,
            ).reshape(natoms, 3)
        except ValueError as exc:
            raise ValueError(f"invalid numeric coordinate record at {path}:{line_number}") from exc
        records.append((iteration, time, positions[projectile_index].copy()))
    if not records:
        raise ValueError(f"no coordinate records found in {path}")
    iterations = np.asarray([item[0] for item in records], dtype=int)
    if np.any(np.diff(iterations) <= 0):
        raise ValueError("coordinate iterations must be strictly increasing")
    return CoordinateSeries(
        iterations=iterations,
        times=np.asarray([item[1] for item in records], dtype=float),
        projectile_positions=np.asarray([item[2] for item in records], dtype=float),
    )


def discover_density_frames(stage_dir: str | Path) -> dict[int, Path]:
    """Return one unambiguous NetCDF density file for each Octopus TD iteration."""

    frames: dict[int, Path] = {}
    for directory in sorted(Path(stage_dir).glob("td.*")):
        match = _TD_DIRECTORY.match(directory.name)
        if not match or not directory.is_dir():
            continue
        exact = directory / "density.ncdf"
        candidates = [exact] if exact.is_file() else sorted(directory.glob("density*.ncdf"))
        if len(candidates) != 1:
            raise ValueError(
                f"{directory} must contain exactly one density NetCDF file; found {len(candidates)}"
            )
        frames[int(match.group(1))] = candidates[0]
    if not frames:
        raise ValueError(f"no td.*/density*.ncdf frames found below {stage_dir}")
    return frames


def integrate_density_on_rectilinear_grid(
    density: np.ndarray,
    axes: tuple[np.ndarray, np.ndarray, np.ndarray],
    center: np.ndarray,
    radius: float,
) -> float:
    """Integrate one density frame in a moving sphere without rescaling."""

    values = np.asarray(density, dtype=float)
    coordinates = tuple(np.asarray(axis, dtype=float) for axis in axes)
    expected = tuple(axis.size for axis in coordinates)
    if values.shape != expected:
        raise ValueError(f"density shape {values.shape} does not match grid shape {expected}")
    spacings = []
    for axis in coordinates:
        if axis.ndim != 1 or axis.size < 2 or np.any(np.diff(axis) <= 0):
            raise ValueError("density axes must be strictly increasing one-dimensional arrays")
        delta = np.diff(axis)
        if not np.allclose(delta, delta[0], rtol=1.0e-10, atol=1.0e-12):
            raise ValueError("nonuniform density grids require explicit quadrature weights")
        spacings.append(float(delta[0]))
    mesh = np.meshgrid(*coordinates, indexing="ij")
    points = np.column_stack([component.reshape(-1) for component in mesh])
    return moving_sphere_population(
        values.reshape(-1), points, np.asarray(center, dtype=float), radius, float(np.prod(spacings))
    )


def load_octopus_density(path: str | Path) -> tuple[np.ndarray, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Load a 3-D Octopus density NetCDF file through xarray.

    The variable and its three coordinate axes must be unambiguous. Spin or
    component-resolved files are rejected instead of being silently combined.
    """

    try:
        import xarray as xr
    except ImportError as exc:
        raise RuntimeError("reading Octopus NetCDF output requires xarray and netCDF4") from exc
    with xr.open_dataset(path, decode_cf=False, mask_and_scale=False) as dataset:
        candidates = []
        for name, variable in dataset.data_vars.items():
            squeezed = variable.squeeze(drop=True)
            if squeezed.ndim == 3 and np.issubdtype(squeezed.dtype, np.number):
                candidates.append((name, squeezed))
        named = [item for item in candidates if "density" in item[0].lower()]
        selected = named if len(named) == 1 else candidates
        if len(selected) != 1:
            names = [item[0] for item in candidates]
            raise ValueError(f"cannot identify one 3-D density variable in {path}; candidates={names}")
        _, variable = selected[0]
        axes = []
        for dimension in variable.dims:
            if dimension not in dataset.coords:
                raise ValueError(f"density dimension {dimension!r} has no coordinate variable")
            axis = np.asarray(dataset.coords[dimension].values, dtype=float).reshape(-1)
            axes.append(axis)
        values = np.asarray(variable.values, dtype=float)
    if np.any(~np.isfinite(values)):
        raise ValueError(f"density contains non-finite values: {path}")
    return values, (axes[0], axes[1], axes[2])


def extract_population_series(
    stage_dir: str | Path,
    *,
    natoms: int,
    projectile_index: int,
    radius: float,
    density_loader: Callable[[str | Path], tuple[np.ndarray, tuple[np.ndarray, np.ndarray, np.ndarray]]] = load_octopus_density,
) -> PopulationSeries:
    stage = Path(stage_dir)
    coordinate_path = stage / "td.general" / "coordinates"
    coordinates = parse_td_coordinates(
        coordinate_path, natoms=natoms, projectile_index=projectile_index
    )
    coordinate_index = {int(value): index for index, value in enumerate(coordinates.iterations)}
    frames = discover_density_frames(stage)
    missing = sorted(set(frames) - set(coordinate_index))
    if missing:
        raise ValueError(f"density iterations lack coordinate records: {missing[:8]}")
    iterations = np.asarray(sorted(frames), dtype=int)
    indices = np.asarray([coordinate_index[int(value)] for value in iterations], dtype=int)
    positions = coordinates.projectile_positions[indices]
    populations = []
    for iteration, center in zip(iterations, positions):
        density, axes = density_loader(frames[int(iteration)])
        populations.append(integrate_density_on_rectilinear_grid(density, axes, center, radius))
    source_files = (coordinate_path, *(frames[int(value)] for value in iterations))
    return PopulationSeries(
        iterations=iterations,
        times=coordinates.times[indices],
        projectile_positions=positions,
        populations=np.asarray(populations, dtype=float),
        source_files=tuple(source_files),
    )


def align_isolated_population(
    production_x: np.ndarray,
    isolated_x: np.ndarray,
    isolated_population: np.ndarray,
) -> np.ndarray:
    """Select recorded populations at matching positions, without resampling.

    The absolute 1e-10 tolerance accommodates coordinate serialization only.
    Missing or ambiguous samples are errors.
    """

    target = np.asarray(production_x, dtype=float)
    source = np.asarray(isolated_x, dtype=float)
    values = np.asarray(isolated_population, dtype=float)
    if source.shape != values.shape or source.ndim != 1:
        raise ValueError("isolated positions and populations must be equal one-dimensional arrays")
    if (source.size == 0 or target.ndim != 1 or np.any(~np.isfinite(source))
            or np.any(~np.isfinite(target)) or np.any(~np.isfinite(values))):
        raise ValueError("reference matching requires finite, nonempty source data")
    order = np.argsort(source)
    source = source[order]
    values = values[order]
    if np.any(np.diff(source) <= 0):
        raise ValueError("isolated projectile x positions must be unique")
    indices = []
    for position in target:
        index = int(np.searchsorted(source, position))
        candidates = [i for i in (index - 1, index, index + 1) if 0 <= i < source.size]
        matches = [i for i in candidates if abs(source[i] - position) <= 1.0e-10]
        if len(matches) != 1:
            raise ValueError(f"production position {position} needs exactly one recorded isolated-reference sample")
        indices.append(matches[0])
    return values[np.asarray(indices, dtype=int)]


def extract_paired_plateau(
    production_stage: str | Path,
    isolated_stage: str | Path,
    *,
    radius: float,
    tail_fraction: float,
    production_natoms: int = 19,
    isolated_natoms: int = 1,
    projectile_index: int = 0,
    density_loader: Callable[[str | Path], tuple[np.ndarray, tuple[np.ndarray, np.ndarray, np.ndarray]]] = load_octopus_density,
) -> tuple[float, dict[str, object]]:
    production = extract_population_series(
        production_stage,
        natoms=production_natoms,
        projectile_index=projectile_index,
        radius=radius,
        density_loader=density_loader,
    )
    isolated = extract_population_series(
        isolated_stage,
        natoms=isolated_natoms,
        projectile_index=projectile_index,
        radius=radius,
        density_loader=density_loader,
    )
    aligned = align_isolated_population(
        production.projectile_positions[:, 0],
        isolated.projectile_positions[:, 0],
        isolated.populations,
    )
    mean_loss_series = paired_projectile_mean_loss(production.populations, aligned)
    statistics = plateau_statistics(mean_loss_series, tail_fraction=tail_fraction)
    mean_loss = float(statistics["mean"])
    detachment_sectors(mean_loss)
    files = (*production.source_files, *isolated.source_files)
    provenance = {
        "production_stage": str(Path(production_stage).resolve()),
        "isolated_stage": str(Path(isolated_stage).resolve()),
        "source_files": [
            {"path": str(path.resolve()), "sha256": sha256_file(path)} for path in files
        ],
        "frame_count": int(production.iterations.size),
        "plateau": statistics,
    }
    return mean_loss, provenance


def load_completed_run(run_dir: str | Path) -> tuple[dict, Path]:
    root = Path(run_dir)
    manifest_path = root / "manifest.json"
    status_path = root / "run_status.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    status = json.loads(status_path.read_text(encoding="utf-8"))
    if status.get("state") != "complete":
        raise ValueError(f"run is not complete: {root}")
    input_path = root / "input" / "inp"
    input_text = input_path.read_text(encoding="utf-8")
    if not re.search(r"^\s*UnitsOutput\s*=\s*eV_Angstrom\s*$", input_text, re.MULTILINE):
        raise ValueError(
            f"{input_path} must declare UnitsOutput=eV_Angstrom for the extraction radius contract"
        )
    stage = root / "stages" / "production_ehrenfest"
    if not stage.is_dir():
        raise FileNotFoundError(stage)
    return manifest, stage


def load_mean_loss_grid(path: str | Path) -> tuple[RTDetachment, dict]:
    """Load a production grid and enforce its direct-TDDFT provenance contract."""

    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if payload.get("source_kind") != "octopus_rt_tddft_density":
        raise ValueError(
            "production yields require a grid extracted directly from Octopus RT-TDDFT density output"
        )
    assert_no_scaling_policy(payload.get("data_policy", {}))
    required = {"heights_bohr", "velocities_au", "mean_loss", "runs"}
    missing = required - set(payload)
    if missing:
        raise ValueError(f"TDDFT grid lacks fields: {sorted(missing)}")
    if not payload["runs"]:
        raise ValueError("TDDFT grid has no run provenance")
    hexadecimal = set("0123456789abcdef")
    for index, run in enumerate(payload["runs"]):
        for field in ("production_manifest_sha256", "isolated_manifest_sha256"):
            digest = str(run.get(field, "")).lower()
            if len(digest) != 64 or set(digest) - hexadecimal:
                raise ValueError(f"TDDFT run {index} has an invalid {field}")
        source_files = run.get("provenance", {}).get("source_files", [])
        if not source_files:
            raise ValueError(f"TDDFT run {index} has no density/coordinate source provenance")
        for source in source_files:
            digest = str(source.get("sha256", "")).lower()
            if len(digest) != 64 or set(digest) - hexadecimal:
                raise ValueError(f"TDDFT run {index} has an invalid source-file digest")
    model = RTDetachment.from_mean_loss_grid(
        np.asarray(payload["heights_bohr"], dtype=float),
        np.asarray(payload["velocities_au"], dtype=float),
        np.asarray(payload["mean_loss"], dtype=float),
    )
    return model, payload
