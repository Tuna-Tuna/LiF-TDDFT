"""Co-moving population and plateau analysis."""

from __future__ import annotations

import numpy as np


def moving_sphere_population(
    density: np.ndarray,
    grid_points: np.ndarray,
    projectile_position: np.ndarray,
    radius: float,
    voxel_volume: float,
) -> float:
    density = np.asarray(density, dtype=float).reshape(-1)
    points = np.asarray(grid_points, dtype=float).reshape(-1, 3)
    center = np.asarray(projectile_position, dtype=float).reshape(3)
    if density.size != points.shape[0]:
        raise ValueError("density and grid point counts differ")
    if radius <= 0 or voxel_volume <= 0:
        raise ValueError("radius and voxel_volume must be positive")
    mask = np.einsum("ij,ij->i", points - center, points - center) <= radius**2
    return float(density[mask].sum() * voxel_volume)


def population_timeseries(
    densities: np.ndarray,
    grid_points: np.ndarray,
    projectile_positions: np.ndarray,
    radii: list[float],
    voxel_volume: float,
) -> dict[float, np.ndarray]:
    rho = np.asarray(densities, dtype=float)
    centers = np.asarray(projectile_positions, dtype=float)
    if rho.shape[0] != centers.shape[0]:
        raise ValueError("one projectile position is required per density frame")
    return {
        radius: np.asarray(
            [moving_sphere_population(frame, grid_points, center, radius, voxel_volume)
             for frame, center in zip(rho, centers)],
            dtype=float,
        )
        for radius in radii
    }


def paired_projectile_mean_loss(
    interacting_population: np.ndarray,
    aligned_isolated_population: np.ndarray,
    *,
    entrance_index: int = 0,
) -> np.ndarray:
    """Return the difference between interacting and isolated population loss."""

    interacting = np.asarray(interacting_population, dtype=float)
    isolated = np.asarray(aligned_isolated_population, dtype=float)
    if interacting.ndim != 1 or interacting.shape != isolated.shape:
        raise ValueError("interacting and aligned isolated populations must match")
    if not 0 <= entrance_index < interacting.size:
        raise IndexError("entrance_index is outside the population series")
    return (interacting[entrance_index] - interacting) - (
        isolated[entrance_index] - isolated
    )


def surface_corrected_local_deficit(
    paired_loss: np.ndarray,
    static_surface_population: np.ndarray,
    *,
    entrance_index: int = 0,
) -> np.ndarray:
    """Apply the moving static-surface correction of the moving static-surface correction."""

    loss = np.asarray(paired_loss, dtype=float)
    surface = np.asarray(static_surface_population, dtype=float)
    if loss.ndim != 1 or loss.shape != surface.shape:
        raise ValueError("loss and static-surface population series must match")
    if not 0 <= entrance_index < loss.size:
        raise IndexError("entrance_index is outside the population series")
    return loss + surface - surface[entrance_index]


def plateau_statistics(values: np.ndarray, tail_fraction: float = 0.2) -> dict[str, float]:
    series = np.asarray(values, dtype=float)
    if series.ndim != 1 or series.size < 5:
        raise ValueError("plateau analysis needs at least five samples")
    count = max(3, int(np.ceil(series.size * tail_fraction)))
    tail = series[-count:]
    x = np.arange(count, dtype=float)
    slope = float(np.polyfit(x, tail, 1)[0])
    return {
        "mean": float(tail.mean()),
        "standard_deviation": float(tail.std(ddof=1)),
        "slope_per_sample": slope,
        "samples": count,
    }
