"""Small deterministic classical integrator for force-field trajectories."""

from __future__ import annotations

from typing import Callable

import numpy as np


def integrate_trajectory(
    initial_position: np.ndarray,
    initial_velocity: np.ndarray,
    mass: float,
    time_step: float,
    steps: int,
    force: Callable[[np.ndarray], np.ndarray],
) -> dict[str, np.ndarray]:
    if mass <= 0 or time_step <= 0 or steps < 1:
        raise ValueError("mass, time step, and steps must be positive")
    positions = np.empty((steps + 1, 3), dtype=float)
    velocities = np.empty_like(positions)
    positions[0] = np.asarray(initial_position, dtype=float)
    velocities[0] = np.asarray(initial_velocity, dtype=float)
    acceleration = np.asarray(force(positions[0]), dtype=float) / mass
    for index in range(steps):
        positions[index + 1] = positions[index] + velocities[index] * time_step + 0.5 * acceleration * time_step**2
        next_acceleration = np.asarray(force(positions[index + 1]), dtype=float) / mass
        velocities[index + 1] = velocities[index] + 0.5 * (acceleration + next_acceleration) * time_step
        acceleration = next_acceleration
    return {"position": positions, "velocity": velocities}


def integrate_grazing_trajectory(
    initial_x_bohr: float,
    initial_height_bohr: float,
    parallel_velocity_au: float,
    normal_velocity_au: float,
    mass_electron_masses: float,
    time_step_au: float,
    steps: int,
    normal_force: Callable[[float, float], float],
) -> dict[str, np.ndarray]:
    """Integrate manuscript Sec. 2.6 trajectories.

    Parallel motion is uniform, while the surface-normal coordinate follows
    velocity Verlet under the velocity-specific force surface ``Fz_v(x,z)``.
    """

    if mass_electron_masses <= 0 or time_step_au <= 0 or steps < 1:
        raise ValueError("mass, time step, and steps must be positive")
    if parallel_velocity_au == 0:
        raise ValueError("parallel velocity must be nonzero")
    time = np.arange(steps + 1, dtype=float) * float(time_step_au)
    x = float(initial_x_bohr) + float(parallel_velocity_au) * time
    z = np.empty(steps + 1, dtype=float)
    vz = np.empty(steps + 1, dtype=float)
    z[0] = float(initial_height_bohr)
    vz[0] = float(normal_velocity_au)
    acceleration = float(normal_force(x[0], z[0])) / float(mass_electron_masses)
    for index in range(steps):
        z[index + 1] = z[index] + vz[index] * time_step_au + 0.5 * acceleration * time_step_au**2
        next_acceleration = float(normal_force(x[index + 1], z[index + 1])) / float(
            mass_electron_masses
        )
        vz[index + 1] = vz[index] + 0.5 * (acceleration + next_acceleration) * time_step_au
        acceleration = next_acceleration
    return {
        "time": time,
        "x": x,
        "height": z,
        "v_parallel": np.full(steps + 1, parallel_velocity_au, dtype=float),
        "v_normal": vz,
    }
