#!/usr/bin/env python3
"""Generate synthetic example data for LiF-TDDFT testing."""
import os
import numpy as np

# Create directory structure
os.makedirs('data/example/h5_0/v0_1', exist_ok=True)

np.random.seed(42)
num_time_steps = 300
num_atoms = 18
num_cols = 1 + 3 * num_atoms + 3 * num_atoms  # time + pos*18 + vel*18 = 109

time_steps = np.arange(1, num_time_steps + 1) * 100

coordinates = np.zeros((num_time_steps, num_cols))
coordinates[:, 0] = time_steps  # time column

# Projectile z position
z_initial = 7.0
z_velocity = -0.1
coordinates[:, 5] = z_initial + z_velocity * time_steps / 100.0

# Small random x, y
coordinates[:, 3] = 0.001 * np.random.randn(num_time_steps)
coordinates[:, 4] = 0.001 * np.random.randn(num_time_steps)

# Velocities for other atoms
for i in range(num_atoms):
    col_start = 1 + 3*num_atoms + 3*i
    coordinates[:, col_start:col_start+3] = 0.01 * np.random.randn(num_time_steps, 3)

# Bound density values
hsemt = 9.0 - 0.2 * np.log(np.arange(1, num_time_steps+1)) + 0.01 * np.random.randn(num_time_steps)
hsemt = np.maximum(hsemt, 5.0)
esemt = 8.5 - 0.3 * np.log(np.arange(1, num_time_steps+1)) + 0.02 * np.random.randn(num_time_steps)
esemt = np.maximum(esemt, 4.0)

# Save as .mat
from scipy.io import savemat
savemat('data/example/h5_0/v0_1/matlab.mat', {
    'coordinates': coordinates,
    'hsemt': hsemt.reshape(1, -1),
    'esemt': esemt.reshape(1, -1),
    'time_steps': time_steps.reshape(-1, 1),
    'bound_density_min': hsemt.reshape(1, -1),
    'bound_density_final': esemt.reshape(1, -1),
})
print(f'Saved matlab.mat: coordinates={coordinates.shape}')
print(f'  hsemt: [{hsemt.min():.2f}, {hsemt.max():.2f}]')
print(f'  esemt: [{esemt.min():.2f}, {esemt.max():.2f}]')

# Write README
readme = """\
# Example Data

This directory contains **synthetic example data** for testing
the LiF-TDDFT analysis pipeline.

**Note:** These data are artificially generated and do NOT
represent real physics results. They are provided solely to
verify that the MATLAB scripts run correctly.

## Structure

```
data/example/
└── h5_0/              # Height = 5.0 au
    └── v0_1/          # Velocity = 0.1 au
        └── matlab.mat # Synthetic TDDFT results
```

## Variables in matlab.mat

| Variable | Shape | Description |
|----------|-------|-------------|
| `coordinates` | [300 x 109] | Trajectory coordinates |
| `hsemt` | [1 x 300] | Bound electron count minimum |
| `esemt` | [1 x 300] | Bound electron count at final time |
| `time_steps` | [300 x 1] | Simulation time steps |
| `bound_density_min` | [1 x 300] | Bound density minimum values |
| `bound_density_final` | [1 x 300] | Bound density final values |

## Usage

```matlab
addpath('matlab/utils');
addpath('matlab/analysis');
load('data/example/h5_0/v0_1/matlab.mat');
```
"""
with open('data/example/README.md', 'w') as f:
    f.write(readme)
print('Saved data/example/README.md')
