# Reproduction workflow

## Supply the actual calculation inputs

Fill `heights_bohr` in `config/campaign_revision.yaml` with recorded calculation
nodes in [1.2,10] bohr. Supply the CAP magnitude for each target velocity in
the manuscript range [0.2,0.4], and actual Troullier-Martins pseudopotential
filenames. The distribution deliberately leaves these fields unset where the
final manuscript does not specify their values. Null fields fail before input
generation. This configuration is not a complete numerical reproduction record.

The central active surface site is at x=0 and the projectile starts at +9
Angstrom. The surface is fixed; ramp segments prescribe each step velocity,
and the production stage maintains the target velocity. The timestep is
0.0004 fs. The grid, box and CAP thickness are 0.1, 16 and 4 Angstrom.

Generate plans with `octopus/generate_inp_files.py`, then execute the independent
paired interacting/isolated histories using `python/run_octopus_workflow.py`.
Dry-run validation does not establish that physical calculations have completed.

## Extract populations and propagate charge states

`scripts/extract_tddft_mean_loss.py` uses the Appendix B smooth weight, with
R=3.5 Angstrom and transition half-width 0.3 Angstrom by default. The same weight
is used for interacting and isolated density frames. The isolated populations
must be recorded at matching positions; reference interpolation remains
disabled as explicitly requested by the author. Coordinates and volume units
must correspond to the configured Angstrom output.

The paired outgoing plateau is Pdet=Ndet. Values outside [0,1] fail instead of
being clipped or converted into a multiple-electron model. Platform means are
direct averages; there is no slope fit or missing-node filling.

`scripts/export_detachment_probabilities.py` exports calculated Pdet nodes
unchanged. Intermediate heights use shape-preserving piecewise cubic PCHIP
at a matching calculated velocity. Extrapolation and velocity interpolation
are rejected. No force, density, capture or final-yield interpolation is added.

Set gamma_c and its source explicitly in `config/demkov_parameters.yaml`.
Supply calculated encounter energy defects in Hartree. The code evaluates Eq.
(9) directly and applies the two-state recursion Eq. (10). It does not estimate
gamma_c from binding energies or infer missing energy-defect values.

The geometry-estimation density script has been removed. The retained
`compare_on_native_grids.m` accepts explicit matching coordinate axes and
densities. There is no automatic atom localization or guessed grid spacing.

The retained MATLAB image-potential interface is
`calc_image_potential(z, v, E0_ev, Ed_ev)` with both dielectric energies required.
Its original integration algorithm, including its known dimension issue,
has not been changed. Do not interpret the interface change as successful
verification of the complete image-potential calculation.

## Validation

Run `python -m unittest discover -s tests -v` and
`python scripts/audit_repository_payload.py`. The tests use explicitly labeled
test inputs and do not calculate manuscript results. For MATLAB/Octave native
grid checks, run `addpath('tests'); test_native_grids`.
`scripts/audit_no_scaling.py` compares exported nodes to their calculated source.

Numerical source datasets and generated outputs remain outside Git. Preserve
their provenance and hashes with the actual results. The code changes alone
do not establish reproduction of the article's figures.
