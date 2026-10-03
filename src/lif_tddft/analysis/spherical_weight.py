"""Smooth spherical weight and analytical gradient from PRA Appendix B."""
import numpy as np


def spherical_weight(points, center, radius, half_width):
    """All coordinates, radius and half_width must use the same length unit."""
    points = np.asarray(points, dtype=float).reshape(-1, 3)
    center = np.asarray(center, dtype=float).reshape(3)
    radius, half_width = float(radius), float(half_width)
    if not np.isfinite(radius) or not np.isfinite(half_width) or not 0 < half_width < radius:
        raise ValueError('smooth sphere requires finite 0 < half_width < radius')
    if np.any(~np.isfinite(points)) or np.any(~np.isfinite(center)):
        raise ValueError('coordinates must be finite')
    relative = points-center
    distance = np.linalg.norm(relative, axis=1)
    w = np.zeros(distance.size)
    gradient = np.zeros_like(relative)
    w[distance <= radius-half_width] = 1.0
    shell = (distance > radius-half_width) & (distance < radius+half_width)
    phase = np.pi*(distance[shell]-radius)/(2*half_width)
    w[shell] = 0.5*(1-np.sin(phase))
    radial_derivative = -np.pi/(4*half_width)*np.cos(phase)
    gradient[shell] = (radial_derivative/distance[shell])[:, None]*relative[shell]
    return w, gradient
