%% =============================================================================
%  fit_effective_potential - Extract equivalent local potential for F⁻ ion from
%                            Octopus DFT ground-state density/potential
%
%  This script reads the static ground-state density and Kohn-Sham potential
%  of a single F⁻ ion from Octopus NetCDF output, performs spherically-averaged
%  3D analysis, and determines the parameter γ (gamma) in the model potential:
%
%    V_model(r) = -14.3996 * erf(γ_Å r) / r    [eV, Å]
%    V_model(r) = -erf(γ r) / r                [Hartree, Bohr]
%
%  THREE METHODS are used to determine γ:
%    Method I:   Density second moment <r²>
%    Method II:  Nonlinear χ² fit of the model potential to v_KS (recommended)
%    Method III: Gaussian tail fit of density (asymptotic region)
%
%  Applicability: Equivalent potential valid for r > 0.8 Å (~1.5 Bohr).
%                 Core region (r < 0.5 Å) must use original pseudopotential.
%
%  OUTPUT:
%    Fminus_local_potential_analysis.png  - 3-panel visualization
%    Console summary with recommended γ parameter
%
%  USAGE:
%    1. Set nc_density_path and nc_potential_path in CONFIGURATION
%    2. Octave users: ensure netcdf package is loaded (see DEPENDENCIES)
%    3. Run the script
%
%  DEPENDENCIES:
%    - MATLAB / GNU Octave 4.4+
%    - Octopus NetCDF density and vks output files
%
%  Octave USERS: ensure netcdf package is loaded:
%    pkg install -forge netcdf
%    pkg load netcdf
%    Or run: run('matlab/utils/octave_startup.m');
%
%  REFERENCES:
%    Project: "Occupied-Space Constraints and Finite-Time Recovery in
%             Electron Detachment during F-/LiF(100) Scattering"
%    Code: Octopus DFT (https://octopus-code.org)
%    Model: V(r) = -erf(γ r)/r  (equivalent local potential approximation)
% =============================================================================

clear; clc; close all;

%% ==================== CONFIGURATION ====================
% Path to F⁻ ground-state density (USER: adjust paths)
nc_density_path = 'Z:\oneF-pbe\static\density.ncdf';

% Path to F⁻ ground-state Kohn-Sham potential (USER: adjust paths)
nc_potential_path = 'Z:\oneF-pbe\static\vks.ncdf';

% Simulation box radius (Å, must match Octopus BoxShape radius)
box_radius = 8;

% Grid spacing (Å, must match Octopus Spacing)
grid_spacing = 0.1;

% Number of radial bins for spherical averaging
num_radial_bins = 400;
% ========================================================

%% ---- 1. Read NetCDF density and potential ----
density = ncread(nc_density_path, 'rdata');
kohn_sham_potential = ncread(nc_potential_path, 'rdata');

% Coordinate grid (Å, consistent with Octopus UnitsOutput = ev_angstrom)
grid_x = -box_radius:grid_spacing:box_radius;
grid_y = -box_radius:grid_spacing:box_radius;
grid_z = -box_radius:grid_spacing:box_radius;

fprintf('Density size: %s\n', mat2str(size(density)));
fprintf('Potential size: %s\n', mat2str(size(kohn_sham_potential)));
fprintf('Grid: %d x %d x %d, Box radius: %.2f Å\n', ...
    length(grid_x), length(grid_y), length(grid_z), max(grid_x));

%% ---- 2. 3D spherical averaging (Å system) ----
[Grid_X, Grid_Y, Grid_Z] = meshgrid(grid_x, grid_y, grid_z);
if ~isequal(size(density), size(Grid_X))
    density = permute(density, [2,1,3]);
    kohn_sham_potential = permute(kohn_sham_potential, [2,1,3]);
end

radial_dist = sqrt(Grid_X.^2 + Grid_Y.^2 + Grid_Z.^2);

radius_max = max(grid_x(:));
radius_edges = linspace(0, radius_max, num_radial_bins + 1);
radius_centers = 0.5 * (radius_edges(1:end-1) + radius_edges(2:end));

density_radial = zeros(num_radial_bins, 1);
potential_radial = zeros(num_radial_bins, 1);
npts = zeros(num_radial_bins, 1);

for i_radial = 1:num_radial_bins
    mask = (radial_dist >= radius_edges(i_radial)) & (radial_dist < radius_edges(i_radial + 1));
    if any(mask(:))
        density_radial(i_radial) = mean(density(mask), 'all');
        potential_radial(i_radial) = mean(kohn_sham_potential(mask), 'all');
        npts(i_radial) = sum(mask(:));
    end
end

valid = npts > 0;
radius_vec   = radius_centers(valid)';
density_vec = density_radial(valid);
potential_vec = potential_radial(valid);

% Force column vectors
radius_vec   = radius_vec(:);
density_vec = density_vec(:);
potential_vec = potential_vec(:);

% Exclude r=0 singularity
if radius_vec(1) == 0
    radius_vec(1) = radius_vec(2) * 0.1;
end

fprintf('Valid radial points: %d, r range: %.2f - %.2f Å\n', ...
    length(radius_vec), min(radius_vec), max(radius_vec));

%% ---- 3. Method I: Density second moment (Å system) ----
dr = mean(diff(radius_vec));
radial_charge = 4 * pi * radius_vec.^2 .* density_vec;   % e/Å

N_total = sum(radial_charge) * dr;
fprintf('\nIntegrated total electron count: %.3f (F⁻ pseudopotential valence ≈ 8.0)\n', N_total);

% Electron count self-check
if abs(N_total - 8.0) > 1.0
    warning('Integrated electron count %.2f deviates from 8.0, check Octopus density units', N_total);
end

% Truncate at 3.2 Å to suppress numerical tail (3.2 Å ≈ 6 Bohr)
cutoff_mask = radius_vec < 3.2;
r2_mean_A = sum(radius_vec(cutoff_mask).^2 .* radial_charge(cutoff_mask)) * dr ...
            / sum(radial_charge(cutoff_mask) * dr);
gamma_A_moment = sqrt(3 / (2 * r2_mean_A));
gamma_moment = gamma_A_moment * 0.529177;   % Convert to Bohr⁻¹

fprintf('Method I (density moment): <r²> = %.4f Å², γ = %.4f Å⁻¹ (%.4f Bohr⁻¹)\n', ...
    r2_mean_A, gamma_A_moment, gamma_moment);

%% ---- 4. Method III: Density tail Gaussian fit (Å system) ----
% Select pure 2p⁶ exponential decay region, excluding core and near-boundary noise
tail_mask = (radius_vec > 0.8) & (radius_vec < 2.5) & (density_vec > 1e-7);

if sum(tail_mask) > 10
    p = polyfit(radius_vec(tail_mask).^2, log(density_vec(tail_mask)), 1);
    gamma_A_tail = sqrt(abs(p(1)));         % Å⁻¹
    gamma_tail = gamma_A_tail * 0.529177;   % Bohr⁻¹
    
    fprintf('Method III (tail Gaussian): γ = %.4f Å⁻¹ (%.4f Bohr⁻¹)\n', ...
        gamma_A_tail, gamma_tail);
else
    gamma_tail = NaN; gamma_A_tail = NaN;
    disp('Method III skipped: insufficient valid tail data');
end

%% ---- 5. Method II: Potential nonlinear fit (Å/eV system) ----
% Model: V_model(r) = -14.3996 * erf(γ_Å r) / r   [eV]
const_eV = 14.3996;

% Hard cutoff core region: erf model cannot describe deep pseudopotential well
weights = ones(size(radius_vec));
weights(radius_vec < 0.8) = 0.0;           % 0.8 Å ≈ 1.5 Bohr
weights(radius_vec > 3.0 & radius_vec < 8.0) = 5.0;
weights(radius_vec > 8.0) = 0.0;

% Scalar objective function (compatible with older MATLAB)
chi2 = @(g) sum(weights .* (potential_vec + const_eV * erf(g * radius_vec) ./ radius_vec).^2, 'all');

% Initialize the optimizer from the calculated density diagnostics.
if ~isnan(gamma_A_tail) && gamma_A_tail > 0.1 && gamma_A_tail < 2.0
    gamma0 = gamma_A_tail;
elseif isfinite(gamma_A_moment) && gamma_A_moment > 0.1 && gamma_A_moment < 2.0
    gamma0 = gamma_A_moment;
else
    error('No valid density-derived starting value for the potential fit.');
end

options = optimset('Display', 'off', 'MaxFunEvals', 1000, 'TolX', 1e-6);
gamma_A_fit = fminsearch(chi2, gamma0, options);
gamma_fit = gamma_A_fit * 0.529177;   % Convert to Bohr⁻¹

% Fit quality (evaluated in eV system)
v_model_eV = -const_eV * erf(gamma_A_fit * radius_vec) ./ radius_vec;
RMSE_eV = sqrt(mean((potential_vec - v_model_eV).^2));
R2 = 1 - sum((potential_vec - v_model_eV).^2) / sum((potential_vec - mean(potential_vec)).^2);

fprintf('Method II (potential fit): γ = %.4f Å⁻¹ (%.4f Bohr⁻¹), RMSE=%.4f eV, R²=%.4f\n', ...
    gamma_A_fit, gamma_fit, RMSE_eV, R2);

%% ---- 6. Validity check and recommended γ ----
if ~isfinite(gamma_fit) || ~isfinite(R2) || gamma_fit > 2.0 || gamma_fit < 0.1 || R2 < 0.85
    error('Potential fit failed its quality criteria; no substitute parameter is selected.');
end
gamma_recommended = gamma_fit;
gamma_A_rec = gamma_A_fit;

%% ---- 7. Visualization ----
figure('Position', [100 100 1500 500]);

% Density comparison (Å system)
subplot(1,3,1);
semilogy(radius_vec, density_vec, 'b-', 'LineWidth', 1.5); hold on;
density_gaussian = N_total * (gamma_A_rec^2 / pi)^(3/2) * exp(-gamma_A_rec^2 * radius_vec.^2);
semilogy(radius_vec, density_gaussian, 'r--', 'LineWidth', 1.5);
xlabel('r (Å)'); ylabel('\rho(r) (e/Å³)');
title('F^- Electron density radial distribution');
legend('Octopus DFT', sprintf('Gaussian (\\gamma=%.3f Å^{-1})', gamma_A_rec), ...
    'Location', 'best');
grid on;

% Potential global comparison (eV)
subplot(1,3,2);
plot(radius_vec, potential_vec, 'b-', 'LineWidth', 2); hold on;
plot(radius_vec, v_model_eV, 'r--', 'LineWidth', 2);
plot(radius_vec, -const_eV ./ radius_vec, 'k:', 'LineWidth', 1.2);
xlabel('r (Å)'); ylabel('V(r) (eV)');
title('F^- Equivalent Local Potential (Global)');
legend('v_{KS}^{DFT}', sprintf('-14.4 erf(%.3fr)/r', gamma_A_rec), ...
    '-14.4/r', 'Location', 'best');
grid on;

% Potential asymptotic region zoom
subplot(1,3,3);
idx = radius_vec > 1.5;
plot(radius_vec(idx), potential_vec(idx), 'b-', 'LineWidth', 2); hold on;
plot(radius_vec(idx), v_model_eV(idx), 'r--', 'LineWidth', 2);
plot(radius_vec(idx), -const_eV ./ radius_vec(idx), 'k:', 'LineWidth', 1.2);
xlabel('r (Å)'); ylabel('V(r) (eV)');
title('Asymptotic region (r > 1.5 Å)');
legend('v_{KS}^{DFT}', sprintf('-14.4 erf(%.3fr)/r', gamma_A_rec), ...
    '-14.4/r', 'Location', 'best');
grid on;

saveas(gcf, 'Fminus_local_potential_analysis.png');
fprintf('\nFigure saved: Fminus_local_potential_analysis.png\n');

%% ---- 8. Results summary ----
fprintf('\n========== F^- Equivalent Local Potential Parameter Summary ==========\n');
fprintf('Analytic form (Å/eV):  V_{F^-}(r) = -14.3996 * erf(\\gamma_A r) / r   [eV, Å]\n');
fprintf('Analytic form (a.u.):   V_{F^-}(r) = -erf(\\gamma r) / r   [Hartree, Bohr]\n');
fprintf('-------------------------------------------\n');
fprintf('\\gamma (density moment):   %.4f Å⁻¹  (%.4f Bohr⁻¹)  [reference only]\n', ...
    gamma_A_moment, gamma_moment);
fprintf('\\gamma (potential fit):    %.4f Å⁻¹  (%.4f Bohr⁻¹)\n', ...
    gamma_A_fit, gamma_fit);
if ~isnan(gamma_tail)
    fprintf('\\gamma (tail Gaussian):    %.4f Å⁻¹  (%.4f Bohr⁻¹)\n', ...
        gamma_A_tail, gamma_tail);
end
fprintf('-------------------------------------------\n');
fprintf('[Recommended] \\gamma = %.4f Bohr⁻¹ (%.4f Å⁻¹)\n', ...
    gamma_recommended, gamma_A_rec);
fprintf('-------------------------------------------\n');
fprintf('Physical check (atomic units):\n');
fprintf('  r \\rightarrow \\infty:  V(r) \\rightarrow -1/r  (correct)\n');
fprintf('  r \\rightarrow 0:      V(0) = -2\\gamma/\\sqrt{\\pi} = %.4f Ha = %.4f eV\n', ...
    -2 * gamma_recommended / sqrt(pi), ...
    -2 * gamma_recommended / sqrt(pi) * 27.2114);
fprintf('  Fit quality: R² = %.4f, RMSE = %.4f eV\n', R2, RMSE_eV);
fprintf('============================================\n');

fprintf('\n[Academic Boundary Statement]\n');
fprintf('This equivalent potential is only valid for r > 0.8 Å (~1.5 Bohr)\n');
fprintf('in the valence and asymptotic regions. The core region (r < 0.5 Å)\n');
fprintf('must retain the original pseudopotential (F.oncvpsp.psp8) local part.\n');
fprintf('Do NOT use -erf(\\gamma r)/r alone to replace the core region.\n');
