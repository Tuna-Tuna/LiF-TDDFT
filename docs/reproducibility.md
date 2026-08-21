# Reproducibility workflow

## 1. Generate RT-TDDFT run plans

`config/campaign_revision.yaml` expands to 8 heights x 6 velocities x 2 paired
variants = 96 independent histories. Each target velocity starts from its own
ground state and segmented acceleration history. Generated `runs/` data are
ignored by Git because full NetCDF outputs are very large.

LiF coordinates use the optional fifth `%Coordinates` column to mark all
surface atoms immobile and the projectile movable. Ground-state stages use
`MoveIons=no`; ramp and production stages use Ehrenfest dynamics with
`MoveIons=yes` and `IonsConstantVelocity=no`.

## 2. Produce the mean-loss map

Run `scripts/extract_tddft_mean_loss.py` on the completed campaign directory.
The production extractor reads each Octopus `td.*/density*.ncdf` frame and the
recorded instantaneous projectile coordinate from `td.general/coordinates`.
It computes projectile-centered populations for the interacting and isolated
projectile runs and applies main Eq. (5) / SI Eq. (S3a) after bounded
position-alignment. The selected population radius must be supplied explicitly
and justified by the radius-convergence calculation.

No out-of-domain node may be silently extrapolated. Local source tables must
label any continuation boundary explicitly. Numerical source tables and their
exported figure values are excluded from Git.

Probability amplitudes are never normalized, rescaled, clipped, or adjusted
by velocity. `scripts/export_figure4_probabilities.py` reads the direct TDDFT
grid, applies only the manuscript sector mapping `Pdet=1-P0`, and writes the
Figure 4 nodes. In particular, `v=0.10` has no special branch. Both the
extracted grid and generated table stay outside Git. In the current manuscript
Figure 4 is `Pdet`, not capture.

## 3. Reconstruct trajectories and encounters

Use `integrate_grazing_trajectory` for uniform parallel motion and calculated
normal motion. `trajectory.events.identify_events` creates a chronological
encounter table. A production table for `compute_sequential_yields.py` must
contain `trajectory_id,event_index,event_time,surface_height,v_parallel,energy_defect_au`.
The detachment mean loss is evaluated from the extracted TDDFT grid inside the
yield program; it must not be supplied as a manually prepared event column.

## 4. Calculate capture and propagate yields

Evaluate SI Eqs. (S15)-(S17), then Eq. (S14). Invoke
`compute_sequential_yields.py` with `--tddft-grid` pointing to the direct
Octopus extraction. Map each mean loss to
`(P0,P1,P2)` and propagate F-/F0/F+ states from neutral F. No experimental
normalization, scale factor, or fit parameter is applied.

Saved MATLAB workspaces and arrays extracted from them are non-production
regression material. They cannot supply the TDDFT grid or the event mean loss.

For the capture panel, `scripts/compute_capture_probability_table.py` evaluates
every row directly with Eq. (S14). The current manuscript labels this quantity
as Figure 7(a) `Pcap`; there is no `v=0.10` multiplier or post-computation
amplitude correction. The loader derives
`gamma=(sqrt(2*Et)+sqrt(2*Ep))/2` from the declared binding-energy inputs,
using `Ep=d_EF` and `Et=d_EF+V_Mad` after conversion to Hartree. It rejects a
directly configured or saved `gamma`/`gama` value.

## 5. Verify

Run the unit suite and `scripts/audit_repository_payload.py` before publishing.
With the required local numerical archives present, also run
`scripts/audit_no_scaling.py` and `scripts/audit_manuscript_consistency.py`.
Generated audit outputs remain ignored by Git. The formula audit is expected
to report a blocker for the archived capture array until the paper and
regenerated Figure 7 source data share one formula.
