# LiF-TDDFT

Code and calculation inputs for the PRA study of electron detachment during
F-/LiF(100) scattering. The repository contains Octopus automation, analysis
of calculated densities and forces, physical capture models, and ordered
charge-state propagation.

## Data handling

- Only detachment probabilities may be interpolated, within a complete
  calculated height/velocity grid. The zero-, one-, and two-electron
  probabilities are obtained at the calculated nodes before interpolation.
- Mean electron loss, force fields, densities, isolated-reference populations,
  capture probabilities, and final yields are not interpolated.
- No data extrapolation, missing-node filling, or substitute numerical
  parameters are used. Missing inputs cause an error.
- Physical models, numerical integration of their equations, and fitting to
  available calculated data are retained. Failed fits do not substitute a
  different estimated result. SciPy is required for special functions.
- Probability amplitudes are not normalized, rescaled, or clipped.

The comparison between full-electronic and point-charge calculations is a
force difference, not by itself a unique decomposition of Pauli energy or
Pauli force. See [model definitions](docs/manuscript_equation_contract.md).

## Installation and calculation workflow

```bash
python -m pip install -e .
python octopus/generate_inp_files.py --config config/campaign_revision.yaml
python python/run_octopus_workflow.py --plan runs/<run_id>/manifest.json --dry-run
```

The dry run checks the plan only. Execute every required paired production
and isolated-projectile run before extraction. Each velocity has an independent
ground-state, acceleration, and production history. Restart failures stop
the run. Supply the actual production pseudopotential files and preserve their
checksums; these files are not distributed here.

```bash
python scripts/extract_tddft_mean_loss.py runs --radius-angstrom <validated-radius> --output results/tddft_mean_loss_grid.json
python scripts/export_detachment_probabilities.py --tddft-grid results/tddft_mean_loss_grid.json
python scripts/compute_sequential_yields.py encounters.csv --tddft-grid results/tddft_mean_loss_grid.json --output results/yields.json
python scripts/compute_capture_probability_table.py <local-input.csv> --output <local-output.csv>
python scripts/audit_no_scaling.py --tddft-grid results/tddft_mean_loss_grid.json
python scripts/audit_repository_payload.py
python -m unittest discover -s tests -v
```

Extraction requires interacting and isolated populations at matching recorded
projectile positions. Differently sampled trajectories must be supplied with
the required matching samples; the code does not resample their populations.
Density comparisons similarly require matching native grids. Force callbacks
must evaluate the supplied physical force law; there is no tabulated-force
interpolation utility.

The mean-loss grid carries direct Octopus density/coordinate provenance.
Saved MATLAB arrays are not accepted as production detachment inputs. The
separate replay utility uses supplied event probabilities only and compares
final yields with experiment only at matching calculated velocities.

See [reproduction instructions](docs/reproducibility.md) and
[local data layout](data/README.md). Numerical datasets, generated figures,
and calculation outputs are outside this source-code repository.

## Software and citation

The workflow targets Python 3.12+, Octopus 16.0, GNU Octave 6.0+ or MATLAB
R2020b+, and the dependencies listed in `pyproject.toml`. MATLAB script input
paths must be configured for the user's calculations. The automated Python
tests do not launch Octopus or establish that production runs are complete.

Use `CITATION.cff` for software citation metadata. License: MIT (`LICENSE`).
