"""Energy-defect components used by the capture model."""

from __future__ import annotations

from collections.abc import Callable

import numpy as np


def _bessel_k0(value: np.ndarray) -> np.ndarray:
    """Evaluate K0 with the required SciPy implementation."""
    from scipy.special import k0

    x = np.asarray(value, dtype=float)
    if np.any(~np.isfinite(x)) or np.any(x <= 0.0):
        raise ValueError("K0 arguments must be finite and positive")
    return np.asarray(k0(x), dtype=float)


def electrostatic_energy_defect(
    projectile_position_bohr: np.ndarray,
    site_positions_bohr: np.ndarray,
    site_charges: np.ndarray,
) -> float:
    """Return the two lattice sums in the energy-defect formula, in Hartree."""

    position = np.asarray(projectile_position_bohr, dtype=float).reshape(3)
    sites = np.asarray(site_positions_bohr, dtype=float)
    charges = np.asarray(site_charges, dtype=float).reshape(-1)
    if sites.shape != (charges.size, 3):
        raise ValueError("site positions and charges are inconsistent")
    reference_r = np.linalg.norm(sites, axis=1)
    projectile_r = np.linalg.norm(position - sites, axis=1)
    if np.any(reference_r == 0.0) or np.any(projectile_r == 0.0):
        raise ValueError("Energy defect is singular at an ionic site")
    return float(np.sum(charges / reference_r) - np.sum(charges / projectile_r))


def mott_littleton_polarization(
    projectile_position_bohr: np.ndarray,
    site_positions_bohr: np.ndarray,
    polarizabilities_bohr3: np.ndarray,
) -> float:
    """Return the induced-dipole energy of the Mott-Littleton polarization formula, in Hartree."""

    position = np.asarray(projectile_position_bohr, dtype=float).reshape(3)
    sites = np.asarray(site_positions_bohr, dtype=float)
    alpha = np.asarray(polarizabilities_bohr3, dtype=float).reshape(-1)
    if sites.shape != (alpha.size, 3):
        raise ValueError("site positions and polarizabilities are inconsistent")
    relative = sites - position
    site_norm = np.linalg.norm(sites, axis=1)
    relative_norm = np.linalg.norm(relative, axis=1)
    if np.any(site_norm == 0.0) or np.any(relative_norm == 0.0):
        raise ValueError("Mott-Littleton polarization is singular at an ionic site")
    # Mott-Littleton polarization uses R_+=0 for the active-site hole and R_-=R for the
    # captured projectile: E_i=(r_i-R_+)/|r_i-R_+|^3
    #                      -(r_i-R_-)/|r_i-R_-|^3.
    field = sites / site_norm[:, None] ** 3 - relative / relative_norm[:, None] ** 3
    return float(-0.5 * np.sum(alpha * np.sum(field**2, axis=1)))


def dynamic_image_interaction(
    v_parallel_au: float,
    height_bohr: float,
    omega_au: np.ndarray,
    surface_response_real: np.ndarray | Callable[[np.ndarray], np.ndarray],
    *,
    projectile_charge: float = 1.0,
) -> float:
    """Numerically evaluate the dynamic-image interaction formula, in Hartree.

    ``projectile_charge`` is the positive charge magnitude ``Q``.
    """

    speed = float(v_parallel_au)
    height = float(height_bohr)
    omega = np.asarray(omega_au, dtype=float)
    charge = float(projectile_charge)
    if speed <= 0.0 or height <= 0.0:
        raise ValueError("v_parallel and Z must be positive")
    if not np.isfinite(charge) or charge <= 0.0:
        raise ValueError("the dynamic-image charge magnitude Q must be positive")
    if omega.ndim != 1 or omega.size < 2 or np.any(np.diff(omega) <= 0.0):
        raise ValueError("omega grid must be a strictly increasing vector")
    response = (
        surface_response_real(omega)
        if callable(surface_response_real)
        else surface_response_real
    )
    response = np.asarray(response, dtype=float)
    if response.shape != omega.shape or np.any(~np.isfinite(response)):
        raise ValueError("surface response must be finite on the omega grid")
    integrand = _bessel_k0(2.0 * omega * height / speed) * response
    trapezoid = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
    return float(-charge / (np.pi * speed) * trapezoid(integrand, omega))


def total_energy_defect(
    electrostatic_hartree: float,
    mott_littleton_hartree: float,
    image_hartree: float,
) -> float:
    """Assemble Delta E from the three contributions in the energy-defect formula."""

    values = np.asarray(
        [electrostatic_hartree, mott_littleton_hartree, image_hartree], dtype=float
    )
    if np.any(~np.isfinite(values)):
        raise ValueError("energy-defect components must be finite")
    return float(values.sum())
