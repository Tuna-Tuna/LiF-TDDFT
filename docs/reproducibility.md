# Reproduction workflow

## Calculate source data

`config/campaign_revision.yaml` defines 8 heights, 6 velocities, and paired
interacting/isolated-projectile variants. Generate the run plans, supply the
actual production pseudopotentials, and complete each ground-state,
acceleration, and Ehrenfest stage. Surface atoms remain fixed; the projectile
is movable. Dry-run plan validation is not a completed calculation.

## Extract recorded populations

`scripts/extract_tddft_mean_loss.py` reads density NetCDF frames and projectile
coordinates directly from completed runs. Supply a population radius justified
by a radius-convergence calculation. The isolated population must have a
recorded sample at every required interacting-projectile position. Matching
uses an absolute coordinate tolerance of 1e-10 only for serialization precision;
it never generates intermediate values. Missing/ambiguous matches fail.

The paired population changes yield a plateau mean loss at each calculated
height/velocity node. Missing grid nodes and nonfinite values fail. Neither
electron counts nor reference populations are interpolated.

## Obtain detachment probabilities

`scripts/export_detachment_probabilities.py` applies the declared sector
mapping to each calculated mean-loss node and exports `Pdet=1-P0` unchanged.
`RTDetachment` separately interpolates the resulting detachment-sector
probabilities inside the complete grid. It rejects all out-of-domain queries.
Mapping the nodes before interpolation is intentional: where a cell spans
the one-electron boundary, this differs from interpolating mean loss and then
applying the piecewise mapping. No mean-loss interpolation is available.

## Evaluate the physical models

Compute the electrostatic, Mott-Littleton, and dynamic-image energy terms at
the requested coordinates. The Demkov formula evaluates capture directly;
it is not interpolated. The capture parameter is derived from binding energies
in `config/demkov_parameters.yaml`. SciPy supplies the Bessel function; absence
of that dependency is an error rather than a reason to use a substitute formula.

Trajectory integration accepts a physical-force callback; no force-grid spline
is supplied. Feed calculated encounter coordinates, velocities, and energy
defects to `scripts/compute_sequential_yields.py` together with the directly
extracted TDDFT grid. Propagation uses chronological F-/F0/F+ transitions.

The density-comparison script requires matching native grids and rejects missing
density values rather than filling them with zero.

## Verify

Run `python -m unittest discover -s tests -v` and
`python scripts/audit_repository_payload.py`. Use `scripts/audit_no_scaling.py`
to compare exported detachment nodes to their calculated sources; optionally
pass `--capture-table` for direct formula verification of a capture table.

Numerical source data and generated outputs remain outside Git. Preserve
input/output hashes and run provenance with local results. The repository
does not synthesize missing production data.

For the native-grid density checks in MATLAB or Octave, run
`addpath('tests'); test_native_grids` from the repository root.
