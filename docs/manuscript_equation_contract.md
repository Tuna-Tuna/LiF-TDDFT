# Definitions aligned with the final PRA manuscript

| Quantity | Implementation | Definition |
| --- | --- | --- |
| Local population | `analysis.spherical_weight`, `analysis.population` | Smooth spherical weight of Appendix B; R=3.5 Angstrom and half-width=0.3 Angstrom for the documented main analysis |
| Reference alignment | `octopus_results.align_isolated_population` | Recorded matching positions only, retained by explicit author instruction; no reference interpolation |
| Detachment | `models.detachment_rt_tddft` | Pdet=Ndet, strictly 0<=Ndet<=1; no clipping |
| Detachment interpolation | `analysis.interpolation` | Shape-preserving piecewise cubic (PCHIP) in height at a recorded velocity; no extrapolation or velocity interpolation |
| Charge propagation | `models.charge_state` | Eq. (10): f_next=(1-f)Pcap+f(1-Pdet), initially neutral |
| Capture | `models.capture_demkov` | Eq. (9), with explicitly supplied gamma_c and energy defect |
| Smooth transport | `analysis.moving_flux.smooth_outward_flux` | -integral[(j-nV) dot grad(w)] dV, with absorption separate |
| Local flyby | `octopus_input`, `generate_inp_files` | Fixed surface, constant prescribed projectile velocity after each ramp increment |
| Grazing trajectory | `trajectory.integrate` | Uniform parallel motion with numerical integration of the normal force |

Appendix B contains the sentence “Reference populations and fluxes are
interpolated along the actual abscissa.” The author explicitly requested
that the existing recorded-sample matching implementation remain unchanged;
this is a known difference from that sentence, not a claim of interpolation.

The manuscript gives height and CAP ranges but not the complete height list,
velocity-dependent CAP assignment, gamma_c value, or actual pseudopotential
filenames. Those inputs are required explicitly. Unknown values are not derived
from unrelated files, estimated or filled with guessed defaults.

The peripheral potential gamma=0.566 is retained because Sec. II A gives it
explicitly. It is distinct from the capture parameter gamma_c.

The MATLAB image-potential algorithm is retained by author instruction. Only
its two dielectric energies are changed to required arguments. This change
does not repair the pre-existing matrix-dimension problem in its integration.
