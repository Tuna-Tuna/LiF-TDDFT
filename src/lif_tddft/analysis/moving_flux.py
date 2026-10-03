"""Moving-control-surface continuity diagnostics."""

from __future__ import annotations

import numpy as np


def outward_flux(
    current_density: np.ndarray,
    electron_density: np.ndarray,
    surface_normals: np.ndarray,
    area_weights: np.ndarray,
    projectile_velocity: np.ndarray,
) -> float:
    current = np.asarray(current_density, dtype=float)
    density = np.asarray(electron_density, dtype=float).reshape(-1)
    normals = np.asarray(surface_normals, dtype=float)
    weights = np.asarray(area_weights, dtype=float).reshape(-1)
    velocity = np.asarray(projectile_velocity, dtype=float).reshape(3)
    if current.shape != normals.shape or current.shape != (density.size, 3):
        raise ValueError("surface current, density, and normal shapes are inconsistent")
    if weights.size != density.size:
        raise ValueError("one area weight is required per surface point")
    relative_current = current - density[:, None] * velocity[None, :]
    return float(np.sum(np.einsum("ij,ij->i", relative_current, normals) * weights))


def closure_residual(times: np.ndarray, populations: np.ndarray, fluxes: np.ndarray, *, absorption_rates: np.ndarray | None = None) -> dict[str, np.ndarray | float]:
    time = np.asarray(times, dtype=float)
    population = np.asarray(populations, dtype=float)
    flux = np.asarray(fluxes, dtype=float)
    if time.ndim != 1 or population.shape != time.shape or flux.shape != time.shape:
        raise ValueError("time, population, and flux must be equal one-dimensional arrays")
    derivative = np.gradient(population, time)
    absorption = np.zeros_like(time) if absorption_rates is None else np.asarray(absorption_rates, dtype=float)
    if absorption.shape != time.shape or np.any(~np.isfinite(absorption)) or np.any(absorption < 0):
        raise ValueError("CAP absorption rates must be finite nonnegative values at every sample")
    residual = derivative + flux + absorption
    scale = max(float(np.max(np.abs(derivative))), float(np.max(np.abs(flux))), 1.0e-15)
    return {
        "dN_dt": derivative,
        "outward_flux": flux,
        "absorption_rate": absorption,
        "residual": residual,
        "relative_rms": float(np.sqrt(np.mean(residual**2)) / scale),
    }


def smooth_outward_flux(current_density, electron_density, grid_points,
                        projectile_position, projectile_velocity,
                        radius, transition_half_width, voxel_volume):
    """Appendix B: -integral[(j-nV) dot grad(w)] dV; CAP loss is separate."""
    from .spherical_weight import spherical_weight
    points = np.asarray(grid_points, dtype=float).reshape(-1, 3)
    density = np.asarray(electron_density, dtype=float).reshape(-1)
    current = np.asarray(current_density, dtype=float)
    velocity = np.asarray(projectile_velocity, dtype=float).reshape(3)
    if current.shape != points.shape or density.size != len(points):
        raise ValueError('density and current must use the same recorded grid')
    if not np.isfinite(voxel_volume) or voxel_volume <= 0:
        raise ValueError('voxel_volume must be positive and finite')
    if np.any(~np.isfinite(density)) or np.any(~np.isfinite(current)) or np.any(~np.isfinite(velocity)):
        raise ValueError('density, current and velocity must be finite')
    _, gradient = spherical_weight(points, projectile_position, radius, transition_half_width)
    return float(-np.sum((current-density[:, None]*velocity)*gradient)*voxel_volume)
