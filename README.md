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
python scripts/compute_sequential_yields.py encounters.csv --output results/yields.json
python scripts/export_figure4_probabilities.py --source <local-table.csv>
python scripts/compute_capture_probability_table.py <local-input.csv> --output <local-output.csv>
python scripts/audit_no_scaling.py
python scripts/audit_manuscript_consistency.py --legacy-dir <local-archive-dir>
python scripts/audit_repository_payload.py
python -m unittest discover -s tests -v
```

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
- `lif_tddft.legacy_run_process`: controlled import and replay of the saved
  `run_process.m` workspace without its network-drive dependencies.

## Importing the saved MATLAB result

`PdE5.mat` is a MATLAB 5 workspace containing the article-result variables.
Extract only the auditable numeric subset and replay its ordered events with:

```bash
matlab -batch "addpath('matlab'); extract_run_process_results('../PdE5.mat','data/processed/run_process')"
python scripts/replay_run_process.py data/processed/run_process
```

The extractor deliberately omits the very large `allCoordinates` cell array,
interpolant objects, temporary structs, and local drive paths.  Python replays
the saved per-event detachment/capture probabilities, reports TD database
coverage, and rejects out-of-domain interpolation rather than silently
extrapolating.

Legacy MATLAB/Octave programs are retained only as regression references.
Production analysis is configuration-driven Python and never reads undeclared
workspace variables or silently fills missing simulation data.

## Non-production P0 evidence analysis

The non-production evidence audit can be rerun without launching Octopus. It
extracts only the compact force-field columns used by the legacy classical
loop, reconstructs ordered collision events and normal velocities, aligns the
continuous Q20 orthogonalization descriptors with the existing v=0.30 TD
trajectory, and compares the archived pre-contract capture array against the
current SI Eq. (S14) implementation:

```bash
matlab -batch "addpath('matlab'); extract_classical_force_fields('../PdE5.mat','data/processed/run_process/force_fields')"
matlab -batch "addpath('matlab'); run_nonproduction_p0_analysis(pwd,'F:/codex/JCTC_纯分析结果与支撑数据_20260727')"
```

Historical reports and their numerical tables remain in the authors' local
archive and are excluded from Git. Current formula-audit outputs can be
produced locally by `scripts/audit_manuscript_consistency.py`. The Wgeom code
remains diagnostic and is not promoted to a formal experimental access weight
until incident-condition, impact-parameter, and trajectory weights are
supplied.

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
