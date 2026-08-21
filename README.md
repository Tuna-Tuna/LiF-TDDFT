# LiF-TDDFT manuscript workflow

This repository supports the manuscript *Occupied-Space Constraints and
Finite-Time Recovery in Electron Detachment during F-/LiF(100) Scattering*
(manuscript ct-2026-013017). The code is scoped to F-/LiF(100) under the
geometries, velocities, cluster models, and electronic-structure settings
declared in the campaign files.

The code does **not** equate the full-electronic minus point-charge difference
with a unique Pauli energy or Pauli force.  Until a frozen-fragment energy
decomposition passes the documented closure gates, outputs use the terms
`fermionic_orthogonalization`, `frozen_fragment_orthogonalization`, and
`force_full_minus_pc`.  The point-charge calculation is an electrostatic
reference without explicit target electrons; it is not a "Pauli-off" model.

The final negative-ion curve is a hybrid result. RT-TDDFT supplies the mean
detached-electron number, the Demkov model supplies capture, and an ordered
charge-state model propagates both along classical trajectories. The mapping
to zero-, one-, and two-electron sectors is explicit; `Pdet = 1-P0`. See
[`docs/manuscript_equation_contract.md`](docs/manuscript_equation_contract.md).

## Release-consistency status

Production code implements SI Eq. (S14) exactly:
`Pcap = 0.5 sech^2[pi(DeltaE+v_parallel^2/2)/(2 gamma v_parallel)]`. It does
not multiply this result by an additional height envelope. The archived local
`zPcapture.csv` predates this contract and contains such a factor. Numerical
archives and generated audit reports are deliberately not tracked by Git.
Authors with the local archive can run
`python scripts/audit_manuscript_consistency.py` to reproduce the discrepancy
report. Public release is blocked until Figure 7 and SI Table S25 are
regenerated with Eq. (S14), or the manuscript explicitly documents a different
capture equation.

The capture parameter is not imported from a saved MATLAB workspace. Production
derives `gamma=(sqrt(2*Et)+sqrt(2*Ep))/2`, with `Ep=d_EF` and
`Et=d_EF+V_Mad` converted to Hartree from the primary inputs declared in
`config/demkov_parameters.yaml`.

Production probability data are never normalized, rescaled, clipped, or
corrected by a velocity-specific factor. This is enforced for Figure 4
detachment probabilities, Figure 7(a) capture probabilities, final yields,
and experimental data. In particular, `v=0.10` has no special amplitude
adjustment. See `config/data_policy.json`.

The exact Li/F pseudopotential files and checksums are also absent from this
working tree. Historical inputs call files named `*.oncvpsp.psp8`, whereas the
manuscript specifies Troullier-Martins pseudopotentials. Do not substitute a
generic library file: public release also requires the authors to identify the
actual production files and reconcile that family label.

## Reproduction levels

1. **example** - generate and validate a small campaign without running
   Octopus: `python octopus/generate_inp_files.py --dry-run`.
2. **local-data reproduction** - analyze author-supplied machine-readable
   tables and rebuild revision figures without committing their numerical
   contents to this repository.
3. **full Octopus reproduction** - execute each target velocity as an
   independent `ground_state -> segmented_velocity_ramp -> production_ehrenfest`
   campaign. Production inputs fix every surface coordinate, allow only the
   projectile to move, and set `IonsConstantVelocity = no`.

## Quick start

```bash
python -m pip install -e .
python octopus/generate_inp_files.py --config config/campaign_revision.yaml
python python/run_octopus_workflow.py --plan runs/<run_id>/manifest.json --dry-run
python scripts/extract_tddft_mean_loss.py runs --radius-angstrom <validated-radius> --output results/tddft_mean_loss_grid.json
python scripts/compute_sequential_yields.py encounters.csv --tddft-grid results/tddft_mean_loss_grid.json --output results/yields.json
python scripts/export_figure4_probabilities.py --tddft-grid results/tddft_mean_loss_grid.json
python scripts/compute_capture_probability_table.py <local-input.csv> --output <local-output.csv>
python scripts/audit_no_scaling.py --tddft-grid results/tddft_mean_loss_grid.json
python scripts/audit_manuscript_consistency.py --legacy-dir <local-archive-dir>
python scripts/audit_repository_payload.py
python -m unittest discover -s tests -v
```

The `--dry-run` command validates a plan only. Remove `--dry-run` and execute
every paired plan before running `extract_tddft_mean_loss.py`; the extractor
requires completed run-status records and refuses incomplete campaigns.

Generated runs are independent by construction.  A target velocity never
reuses the acceleration history of a different target velocity.  Failed or
incomplete restart stages abort the run and are recorded in `run_status.json`.

## Main modules

- `lif_tddft.campaign`: campaign matrix and unique run identifiers.
- `lif_tddft.octopus_input`: input rendering and physical preflight checks.
- `lif_tddft.provenance`: input, pseudopotential, software, Git, and host
  provenance.
- `lif_tddft.analysis`: co-moving populations, moving-control-surface flux,
  orthogonalization descriptors, forces, work/energy, convergence, and bounded
  interpolation.
- `lif_tddft.clusters`: C9/C25 one-active-site topology construction and
  comparison gates.
- `lif_tddft.trajectory`: integration, ordered collision events, and purely
  geometric `W_geom` diagnostics.
- `lif_tddft.models`: mean-loss sector closure, SI Eqs. (S14)-(S17), and
  ordered F-/F0/F+ propagation with the stated F+ -> F0 closure.
- `lif_tddft.octopus_results`: strict extraction of co-moving projectile
  populations and paired mean loss directly from Octopus density NetCDF and
  `td.general/coordinates` output.

## Production TDDFT source contract

Production detachment data are generated only from completed, paired
`production` and `isolated_projectile` Octopus runs. The extractor reads the
time-resolved `td.*/density*.ncdf` files and the instantaneous projectile
positions in `td.general/coordinates`, integrates the same moving sphere in
both runs, position-aligns the isolated reference without extrapolation, and
applies main Eq. (5) / SI Eq. (S3a). It then reports the unscaled plateau mean
on the complete height/velocity grid.

```bash
python scripts/extract_tddft_mean_loss.py runs \
  --radius-angstrom <validated-radius> \
  --output results/tddft_mean_loss_grid.json
```

`compute_sequential_yields.py` requires this extracted grid and rejects files
whose `source_kind` is not `octopus_rt_tddft_density`. It cannot use a saved
MATLAB workspace, archived `zPloss.csv`, or manually populated `mean_loss`
column as production input. Historical MATLAB/Octave replay utilities are
retained only as non-production regression aids and are not part of the
reproduction commands above.

## Data and provenance

Raw NetCDF outputs, experimental values, processed manuscript nodes, extracted
MATLAB arrays, figures, and generated audit tables are intentionally excluded
from Git. They remain local or belong in a separately governed data archive.
Every locally generated table should point to a run manifest containing
SHA-256 hashes of inputs and pseudopotentials, the repository commit, Octopus
version, MPI command, host, and stage status. Missing data are errors;
out-of-domain trajectory events are reported and are not extrapolated.

The public Git payload is limited to source code, tests, calculation inputs,
templates, configuration, and supporting documentation. See `.gitignore` and
run `python scripts/audit_repository_payload.py` before publishing.

## Software

The workflow targets Python 3.12+, Octopus 16.0, GNU Octave 6.0+ or MATLAB
R2020b+, and an MPI implementation appropriate to the target cluster.  Tests
that do not launch Octopus use only Python and NumPy.

## Manuscript status

The manuscript is under revision.  No DOI is assigned here and no publication
claim is implied.  Use `CITATION.cff` for the current author and title metadata.

## License

MIT.  See `LICENSE`.
