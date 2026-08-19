"""Work-energy consistency checks along prescribed coordinates."""

import numpy as np


def cumulative_work(positions: np.ndarray, forces: np.ndarray) -> np.ndarray:
    r = np.asarray(positions, dtype=float)
    f = np.asarray(forces, dtype=float)
    if r.shape != f.shape or r.ndim != 2 or r.shape[1] != 3:
        raise ValueError("positions and forces must have shape (n, 3)")
    increments = np.einsum("ij,ij->i", 0.5 * (f[1:] + f[:-1]), r[1:] - r[:-1])
    return np.concatenate(([0.0], np.cumsum(increments)))


def work_energy_residual(work: np.ndarray, kinetic_energy: np.ndarray) -> np.ndarray:
    work = np.asarray(work, dtype=float)
    kinetic = np.asarray(kinetic_energy, dtype=float)
    if work.shape != kinetic.shape:
        raise ValueError("work and kinetic energy arrays must match")
    return kinetic - kinetic[0] - work
