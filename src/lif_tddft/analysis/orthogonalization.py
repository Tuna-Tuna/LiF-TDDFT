"""Frozen-fragment occupied-space orthogonalization descriptors."""

from __future__ import annotations

import numpy as np


def cross_overlap(orbitals_a: np.ndarray, orbitals_b: np.ndarray, weights: np.ndarray) -> np.ndarray:
    a = np.asarray(orbitals_a, dtype=complex)
    b = np.asarray(orbitals_b, dtype=complex)
    w = np.asarray(weights, dtype=float).reshape(-1)
    if a.ndim != 2 or b.ndim != 2 or a.shape[1] != b.shape[1] or a.shape[1] != w.size:
        raise ValueError("orbitals must be (n_orbitals, n_grid) on one common grid")
    return (a.conj() * w[None, :]) @ b.T


def occupied_overlap_descriptors(s_ab: np.ndarray) -> dict[str, float]:
    singular = np.linalg.svd(np.asarray(s_ab), compute_uv=False)
    omega = float(np.sum(singular**2))
    return {
        "maximum_singular_value": float(singular.max(initial=0.0)),
        "Omega_occ": omega,
    }


def lowdin_orthogonalize(orbitals: np.ndarray, weights: np.ndarray, cutoff: float = 1.0e-10):
    phi = np.asarray(orbitals, dtype=complex)
    w = np.asarray(weights, dtype=float).reshape(-1)
    overlap = (phi.conj() * w[None, :]) @ phi.T
    eigenvalues, eigenvectors = np.linalg.eigh(overlap)
    if eigenvalues.min() <= cutoff:
        raise ValueError("overlap matrix is ill-conditioned at the requested cutoff")
    inverse_sqrt = eigenvectors @ np.diag(eigenvalues ** -0.5) @ eigenvectors.conj().T
    ortho = inverse_sqrt @ phi
    diagnostics = {
        "minimum_eigenvalue": float(eigenvalues.min()),
        "condition_number": float(eigenvalues.max() / eigenvalues.min()),
        "cutoff": cutoff,
    }
    return ortho, diagnostics


def canonical_orthogonalize(orbitals: np.ndarray, weights: np.ndarray, cutoff: float = 1.0e-10):
    phi = np.asarray(orbitals, dtype=complex)
    w = np.asarray(weights, dtype=float).reshape(-1)
    weighted = phi * np.sqrt(w)[None, :]
    _, singular, vh = np.linalg.svd(weighted, full_matrices=False)
    if singular.min() ** 2 <= cutoff:
        raise ValueError("canonical orthogonalization cutoff removes an occupied orbital")
    ortho = vh / np.sqrt(w)[None, :]
    return ortho, {"minimum_singular_value": float(singular.min()), "cutoff": cutoff}


def density_rearrangement(original: np.ndarray, orthogonalized: np.ndarray, occupations: np.ndarray) -> np.ndarray:
    before = np.asarray(original, dtype=complex)
    after = np.asarray(orthogonalized, dtype=complex)
    occ = np.asarray(occupations, dtype=float).reshape(-1, 1)
    if before.shape != after.shape or before.shape[0] != occ.size:
        raise ValueError("orbital and occupation shapes are inconsistent")
    return np.sum(occ * (np.abs(after) ** 2 - np.abs(before) ** 2), axis=0).real
