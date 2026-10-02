# Data Availability

This directory documents the expected local data layout. The Git repository
does not distribute experimental values, simulation outputs, processed tables,
figures, or extracted workspaces. Only source code and calculation inputs are
published; numerical datasets remain local or belong in a separately governed
data archive.

## Data Description

### Octopus Output Structure

A typical run directory contains:

```
run_directory/
├── inp                                        # Octopus input file
├── out                                         # Main Octopus log
├── restart/                                    # Wavefunction restart files
├── static/                                     # Ground-state output
│   └── density.ncdf                            # Ground-state electron density
├── td.general/                                 # Time-dependent general output
│   ├── coordinates                             # Time-resolved atomic coordinates
│   ├── energy                                  # Time-resolved total energy
│   ├── eigenvalues                             # Time-resolved Kohn-Sham eigenvalues
│   └── ion_ch                                  # Ionization channel probabilities
├── td.0000100/                                 # Output interval (every 100 steps)
│   └── density.ncdf                            # Time-resolved 3D electron density
├── td.0000200/
│   └── density.ncdf
├── ...
└── td.000NN00/
    └── density.ncdf
```

### Key Data Files

| File | Format | Description | Typical Size |
|------|--------|-------------|-------------|
| `static/density.ncdf` | NetCDF | Ground-state 3D electron density | ~500 MB |
| `td.general/coordinates` | ASCII text | Projectile trajectory coordinates | ~1-10 MB |
| `td.general/energy` | ASCII text | Time-resolved total energy | ~1 MB |
| `td.general/eigenvalues` | ASCII text | Kohn-Sham eigenvalues over time | ~1-10 MB |
| `td.NNNNNNN/density.ncdf` | NetCDF | Time-resolved 3D density snapshots | ~500 MB each |

### Data Access

Numerical data are not distributed in this code repository. Contact the
manuscript's corresponding authors for access to the underlying calculation
outputs and experimental comparison data.

### Dataset Organization

The simulation was performed for the following parameter combinations:

**Surface heights (h, au):** 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 7.0, 10.0  
**Projectile velocities (v, au):** 0.10, 0.15, 0.20, 0.30, 0.40, 0.50  

Data is organized as:
```
data/
├── h{h_value}/
│   └── v{v_value}/
│       ├── static/
│       ├── td.general/
│       └── td.{step}/density.ncdf
```

## Using the Analysis Scripts

1. Create the expected directory structure with your Octopus output
2. Update the `CONFIGURATION` section in each MATLAB script to point to your data paths
3. Ensure MATLAB utility functions are on the path:
   ```matlab
   addpath('matlab/utils');
   addpath('matlab/analysis');
   ```

## External Dependencies

- **Pseudopotentials**: the inputs reference `Li.oncvpsp.psp8` and
  `F.oncvpsp.psp8`. The actual production files and their checksums must be
  supplied and verified before calculation; they are not distributed here.
- **Octopus code**: Version 16.0 or later from https://octopus-code.org

## Data contract

Raw Octopus output remains outside Git. Each production run should preserve its
input, stage logs, restart validation, `run_status.json`, and `provenance.json`.
Processed population, moving-flux, force, orthogonalization, event,
`W_geom`, and hybrid-yield tables may be written locally under
`data/processed`, which is ignored by Git, or deposited in a separately
governed data archive. Current-manuscript nodes must remain separate from
legacy MATLAB arrays. Missing production data must not be silently
interpolated or represented as calculated results.
