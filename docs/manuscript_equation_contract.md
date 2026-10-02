# Model definitions

The definitions below identify implemented physical quantities without relying
on manuscript-specific figure or equation numbering. A change in the physical
model requires an explicit scientific justification.

| Quantity | Implementation | Definition |
| --- | --- | --- |
| Paired mean loss | `octopus_results.extract_paired_plateau`, `analysis.population.paired_projectile_mean_loss` | Interacting minus isolated population change at matching recorded positions |
| Surface-corrected deficit | `analysis.population.surface_corrected_local_deficit` | Adds the time-dependent moving static-surface correction |
| Detachment sectors | `models.detachment_rt_tddft.detachment_sectors` | For `N<=1`, `(P0,P1,P2)=(1-N,N,0)`; for `1<N<2`, `(0,2-N,N-1)` |
| Detachment interpolation | `models.detachment_rt_tddft.RTDetachment` | Map each calculated mean-loss node to sectors first; interpolate only those probabilities in-domain |
| Charge-state propagation | `models.charge_state` | Chronological independent encounters; transient F+ returns to F0 at the next encounter |
| Demkov capture | `models.capture_demkov.DemkovCapture` | `0.5*sech^2[pi*(DeltaE+v^2/2)/(2*gamma*v)]`, evaluated directly |
| Capture parameter | `models.demkov_parameters` | `gamma=(sqrt(2*Et)+sqrt(2*Ep))/2`, `Ep=d_EF`, `Et=d_EF+V_Mad`, energies in Hartree |
| Energy defect | `models.capture_energy` | Electrostatic lattice sums plus Mott-Littleton and dynamic-image terms |
| Dynamic-image interaction | `models.capture_energy.dynamic_image_interaction` | `-Q/(pi*v)` times the K0-weighted real surface-response integral |
| Trajectory | `trajectory.integrate.integrate_grazing_trajectory` | Uniform parallel motion and numerical integration of normal motion under the supplied physical force |

No probability scaling, clipping, experimental normalization, or additional
capture height envelope is applied. Missing calculated data cause errors;
neither out-of-domain values nor substitute estimates are generated.

All model inputs use atomic units unless a field name specifies otherwise.
Heights are bohr, velocities atomic units, energies Hartree, `gamma` inverse
bohr, and polarizabilities bohr cubed. Statistical means are computed from the
supplied data and are never used to fill missing production measurements or
calculations. Model parameters and their sources are declared in the configuration.
