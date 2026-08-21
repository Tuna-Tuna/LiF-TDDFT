"""Energy-defect components used by SI Eqs. (S15)--(S17)."""

from __future__ import annotations

from collections.abc import Callable

import numpy as np


def _bessel_i0(value: np.ndarray) -> np.ndarray:
    """Vectorized I0 approximation used only when SciPy is unavailable."""

    x = np.abs(np.asarray(value, dtype=float))
    result = np.empty_like(x)
    small = x < 3.75
    y = (x[small] / 3.75) ** 2
    result[small] = 1.0 + y * (
        3.5156229 + y * (3.0899424 + y * (1.2067492 + y * (
            0.2659732 + y * (0.0360768 + y * 0.0045813)
        )))
    )
    y = 3.75 / x[~small]
    result[~small] = np.exp(x[~small]) / np.sqrt(x[~small]) * (
        0.39894228 + y * (0.01328592 + y * (0.00225319 + y * (
            -0.00157565 + y * (0.00916281 + y * (-0.02057706 + y * (
                0.02635537 + y * (-0.01647633 + y * 0.00392377)
            )))
        )))
    )
    return result


def _bessel_k0(value: np.ndarray) -> np.ndarray:
    """Return K0 using SciPy or the standard polynomial fallback."""

    x = np.asarray(value, dtype=float)
    if np.any(x <= 0.0):
        raise ValueError("K0 arguments must be positive")
    try:
        from scipy.special import k0 as scipy_k0

        return np.asarray(scipy_k0(x), dtype=float)
    except ImportError:
        result = np.empty_like(x)
        small = x <= 2.0
        y = x[small] ** 2 / 4.0
        result[small] = -np.log(x[small] / 2.0) * _bessel_i0(x[small]) + (
            -0.57721566 + y * (0.42278420 + y * (0.23069756 + y * (
                0.03488590 + y * (0.00262698 + y * (0.00010750 + y * 0.00000740))
            )))
        )
        y = 2.0 / x[~small]
        result[~small] = np.exp(-x[~small]) / np.sqrt(x[~small]) * (
            1.25331414 + y * (-0.07832358 + y * (0.02189568 + y * (
                -0.01062446 + y * (0.00587872 + y * (-0.00251540 + y * 0.00053208))
            )))
        )
        return result


def electrostatic_energy_defect(
    projectile_position_bohr: np.ndarray,
    site_positions_bohr: np.ndarray,
    site_charges: np.ndarray,
) -> float:
    """Return the two lattice sums in SI Eq. (S15), in Hartree."""

    position = np.asarray(projectile_position_bohr, dtype=float).reshape(3)
    sites = np.asarray(site_positions_bohr, dtype=float)
    charges = np.asarray(site_charges, dtype=float).reshape(-1)
    if sites.shape != (charges.size, 3):
        raise ValueError("site positions and charges are inconsistent")
    reference_r = np.linalg.norm(sites, axis=1)
    projectile_r = np.linalg.norm(position - sites, axis=1)
    if np.any(reference_r == 0.0) or np.any(projectile_r == 0.0):
        raise ValueError("Eq. (S15) is singular at an ionic site")
    return float(np.sum(charges / reference_r) - np.sum(charges / projectile_r))


def mott_littleton_polarization(
    projectile_position_bohr: np.ndarray,
    site_positions_bohr: np.ndarray,
    polarizabilities_bohr3: np.ndarray,
) -> float:
    """Return the induced-dipole energy of SI Eq. (S16), in Hartree."""

    position = np.asarray(projectile_position_bohr, dtype=float).reshape(3)
    sites = np.asarray(site_positions_bohr, dtype=float)
    alpha = np.asarray(polarizabilities_bohr3, dtype=float).reshape(-1)
    if sites.shape != (alpha.size, 3):
        raise ValueError("site positions and polarizabilities are inconsistent")
    relative = sites - position
    site_norm = np.linalg.norm(sites, axis=1)
    relative_norm = np.linalg.norm(relative, axis=1)
    if np.any(site_norm == 0.0) or np.any(relative_norm == 0.0):
        raise ValueError("Eq. (S16) is singular at an ionic site")
    # Eq. (S16) uses R_+=0 for the active-site hole and R_-=R for the
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
    """Numerically evaluate SI Eq. (S17), in Hartree.

    ``projectile_charge`` is the positive magnitude ``Q`` printed in Eq. (S17).
    """

    speed = float(v_parallel_au)
    height = float(height_bohr)
    omega = np.asarray(omega_au, dtype=float)
    charge = float(projectile_charge)
    if speed <= 0.0 or height <= 0.0:
        raise ValueError("v_parallel and Z must be positive")
    if not np.isfinite(charge) or charge <= 0.0:
        raise ValueError("the Eq. (S17) charge magnitude Q must be positive")
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
    """Assemble Delta E from the three contributions in SI Eq. (S15)."""

    values = np.asarray(
        [electrostatic_hartree, mott_littleton_hartree, image_hartree], dtype=float
    )
    if np.any(~np.isfinite(values)):
        raise ValueError("energy-defect components must be finite")
    return float(values.sum())
