# Manuscript equation contract

This file is the code-to-paper contract for manuscript `ct-2026-013017`.
Production functions may change only together with the corresponding main-text
or Supporting Information equation.

| Paper definition | Code | Enforced interpretation |
|---|---|---|
| Main Eq. (5); SI Eq. (S3a) | `octopus_results.extract_paired_plateau` and `analysis.population.paired_projectile_mean_loss` | Direct Octopus density extraction followed by paired interacting-minus-position-aligned-isolated population change |
| SI Eq. (S3b) | `analysis.population.surface_corrected_local_deficit` | Moving static-surface correction is time dependent |
| Main Sec. 2.6 piecewise closure | `models.detachment_rt_tddft.detachment_sectors` | `Nbar_det<=1`: `(P0,P1,P2)=(1-Nbar,Nbar,0)`; `1<Nbar<2`: `(0,2-Nbar,Nbar-1)` |
| Main Eq. (20) | `models.charge_state` | Chronological independent encounters; the one-electron branch reduces to Eq. (20) |
| SI low-height closure | `models.charge_state.transition_matrix` | Two-electron loss creates transient F+; the next encounter imposes `P(F+->F0)=1` |
| SI Eq. (S14) | `models.capture_demkov.DemkovCapture` | `0.5*sech^2[pi*(DeltaE+v^2/2)/(2*gamma*v)]`; no extra height factor |
| SI Eq. (S15) | `models.capture_energy.electrostatic_energy_defect` and `total_energy_defect` | Two lattice sums plus Mott-Littleton and dynamic-image terms |
| SI Eq. (S16) | `models.capture_energy.mott_littleton_polarization` | `-alpha_pm/2` times the squared field difference between active-site hole `R_+=0` and projectile `R_-=R` |
| SI Eq. (S17) | `models.capture_energy.dynamic_image_interaction` | `-Q/(pi*v_parallel)` times the K0-weighted real surface-response integral |
| Main Sec. 2.6 trajectory | `trajectory.integrate.integrate_grazing_trajectory` | Uniform parallel motion; normal motion integrated under `Fz_v(x,z)` |

## No-scaling contract

All production probability tables use the identity data transform:
`normalization=none`, `scale_factor=1.0`, `value_transform=none`, and
`velocity_specific_adjustments=none`. Out-of-range values cause errors; they
are never clipped. This applies to Figure 4 `Pdet`, Figure 7(a) `Pcap`, final
F- fractions, and the experimental comparison. The `v=0.10` series is not a
special case.

## Units

All model inputs are atomic units unless a field name says otherwise. Heights
and positions are bohr, velocities are atomic units, energy defects are
Hartree, `gamma` is inverse bohr, and ionic polarizabilities are bohr cubed.

## Archived-array boundary

The authors' local `data/processed/run_process` directory is a numerical
extraction of the legacy `PdE5.mat` workspace and is excluded from Git. Its
local `zPcapture.csv` is not an evaluation of current SI Eq. (S14): the legacy
expression multiplied Eq. (S14) by
`sech^2[1.2*v^2*(h-3)]`. Production Python forbids that factor. The difference
can be generated locally by `scripts/audit_manuscript_consistency.py`.

No value extracted from `PdE5.mat`, `run_process.m`, `zPloss.csv`, or another
saved result array is a production detachment input. Production mean loss is
recomputed from paired Octopus `td.*/density*.ncdf` and
`td.general/coordinates` outputs by `scripts/extract_tddft_mean_loss.py`.
`scripts/compute_sequential_yields.py` rejects any grid that does not declare
`source_kind=octopus_rt_tddft_density` and carry run provenance.
