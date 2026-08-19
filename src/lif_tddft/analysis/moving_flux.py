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


def closure_residual(times: np.ndarray, populations: np.ndarray, fluxes: np.ndarray) -> dict[str, np.ndarray | float]:
    time = np.asarray(times, dtype=float)
    population = np.asarray(populations, dtype=float)
    flux = np.asarray(fluxes, dtype=float)
    if time.ndim != 1 or population.shape != time.shape or flux.shape != time.shape:
        raise ValueError("time, population, and flux must be equal one-dimensional arrays")
    derivative = np.gradient(population, time)
    residual = derivative + flux
    scale = max(float(np.max(np.abs(derivative))), float(np.max(np.abs(flux))), 1.0e-15)
    return {
        "dN_dt": derivative,
        "outward_flux": flux,
        "residual": residual,
        "relative_rms": float(np.sqrt(np.mean(residual**2)) / scale),
    }
