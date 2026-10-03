# LiF-TDDFT

Code and calculation inputs for the PRA study of electron detachment during
F-/LiF(100) scattering. The repository contains Octopus automation, analysis
of calculated densities and forces, physical capture models, and ordered
charge-state propagation.

## Data handling

- Only detachment probabilities may be interpolated: shape-preserving PCHIP
  in height at a calculated velocity. Pdet equals the corrected outgoing
  population loss in the single-electron approximation, with 0<=Pdet<=1.
- Mean electron loss, force fields, densities, isolated-reference populations,
  capture probabilities, and final yields are not interpolated.
- No data extrapolation, missing-node filling, or substitute numerical
  parameters are used. Missing inputs cause an error.
- Physical equations are evaluated with the declared model parameters.
  SciPy is required for special functions.
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

First fill the required height nodes, per-velocity CAP magnitudes and actual
pseudopotential filenames in the campaign configuration. These values are not
guessed from the ranges in the manuscript. Local production uses a prescribed
constant velocity. The dry run checks the plan only. Execute every required paired production
and isolated-projectile run before extraction. Each velocity has an independent
ground-state, acceleration, and production history. Restart failures stop
the run. Supply the actual production pseudopotential files and preserve their
checksums; these files are not distributed here.

```bash
python scripts/extract_tddft_mean_loss.py runs --radius-angstrom 3.5 --transition-half-width-angstrom 0.3 --output results/tddft_mean_loss_grid.json
python scripts/export_detachment_probabilities.py --tddft-grid results/tddft_mean_loss_grid.json
python scripts/compute_sequential_yields.py encounters.csv --tddft-grid results/tddft_mean_loss_grid.json --output results/yields.json
python scripts/compute_capture_probability_table.py <local-input.csv> --output <local-output.csv>
python scripts/audit_no_scaling.py --tddft-grid results/tddft_mean_loss_grid.json
python scripts/audit_repository_payload.py
python -m unittest discover -s tests -v
```

Set the actual capture parameter and its source in `config/demkov_parameters.yaml`
before capture or yield calculations. The manuscript does not give its numerical
value, so the configuration contains no guessed default.

Extraction uses the Appendix B smooth spherical weight and requires interacting
and isolated populations at matching recorded
projectile positions. Differently sampled trajectories must be supplied with
the required matching samples; the code does not resample their populations.
The automatic geometry-estimation script has been removed. Density comparisons
require explicitly supplied matching native grids. Force callbacks
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

The legacy MATLAB image-potential algorithm is retained with required dielectric
energy arguments; its pre-existing dimension issue is not fixed by this change.

Use `CITATION.cff` for software citation metadata. License: MIT (`LICENSE`).
