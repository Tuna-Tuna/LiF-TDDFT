%% =============================================================================
%  compare_cluster_slab_density - Compare F⁻ electron density between LiF cluster
%                                 and slab environments (with alignment fix)
%
%  This script compares the electron density distribution of a surface F⁻ ion
%  in two computational environments:
%    (a) LiF cluster (finite system, spherical box)
%    (b) LiF slab   (periodic surface, parallelepiped box)
%
%  The comparison reveals the effect of the extended surface environment on the
%  F- electron-cloud polarization and occupied-space deformation.
%
%  Key analysis steps:
%    1. Import density from Octopus NetCDF (convert e/Å³ → e/Bohr³)
%    2. Locate surface F⁻ ion in each environment (auto-detection)
%    3. Align z-coordinates for fair comparison
%    4. Compare recorded, matching density grids centered at F⁻
%    5. Spherical cutoff at R_cut = 2.8 au
%    6. Compute 2D xy-slices, radial profile Δρ(r), z-axis profile Δρ(z)
%    7. Compute charge center of mass (density-weighted centroid)
%
%  OUTPUT FILES:
%    cluster_xy_z0.txt          - (a) Cluster 2D density at z=0
%    slab_xy_z0.txt             - (b) Slab 2D density at z=0
%    radial_density_diff.txt    - (c) Δρ(r) radial profile
%    z_axis_density_diff.txt    - (d) Δρ(z) along surface normal
%    z_axis_cumulative.txt      - Cumulative charge difference vs z
%    charge_center_summary.txt  - Charge centroid and alignment info
%
%  USAGE:
%    1. Set CONFIGURATION paths/parameters
%    2. Octave users: ensure netcdf package is loaded:
%         pkg install -forge netcdf
%         pkg load netcdf
%       Or run: run('matlab/utils/octave_startup.m');
%    3. Run the script
%
%  DEPENDENCIES:
%    - MATLAB / GNU Octave 4.4+
%    - Octopus static/density.ncdf NetCDF files
%
%  Octave USERS: ensure netcdf package is loaded (see octave_startup.m).
%
%  REFERENCES:
%    Project: "Occupied-Space Constraints and Finite-Time Recovery in
%             Electron Detachment during F-/LiF(100) Scattering"
%    Code: Octopus DFT (https://octopus-code.org)
% =============================================================================
clear; clc; close all;

%% ==================== CONFIGURATION ====================
% Physical constants
bohr2ang = 0.5291772083;

% Unit conversion factor: e·Å^{-3} -> e·Bohr^{-3}
convert_Ang3_to_bohr3 = bohr2ang^3;   % = 0.148184

% NetCDF density file paths (USER: adjust to your local paths)
cluster_path = "Z:\version1.5\LiF\v1.4\static\density.ncdf";
slab_path    = "Z:\version1.5\LiF\slab\static\density.ncdf";

% Slab grid spacing (au, USER: adjust to match Octopus Spacing for slab)
slab_spacing_au = 0.20;

% Analysis parameters
cutoff_radius = 2.8;     % au, spherical cutoff for density comparison
box_half_width = 3.5;    % au, local extraction box half-width

% Crystal geometry (LiF rock-salt, surface orientation)
lattice_constant = 3.7958;  % au
% ========================================================

%% ---- 1. Import density data (e·Å^{-3} → e·Bohr^{-3}) ----
density_cluster = ncread(cluster_path, "rdata") * convert_Ang3_to_bohr3;  % [z,y,x]
density_slab    = ncread(slab_path,    "rdata") * convert_Ang3_to_bohr3;  % [z,y,x]

[Nz_c, Ny_c, Nx_c] = size(density_cluster);
[Nz_s, Ny_s, Nx_s] = size(density_slab);

fprintf('=== Data Import ===\n');
fprintf('Cluster dimensions [z,y,x]: [%d, %d, %d]\n', Nz_c, Ny_c, Nx_c);
fprintf('Slab    dimensions [z,y,x]: [%d, %d, %d]\n', Nz_s, Ny_s, Nx_s);
fprintf('Unit conversion: %.6f (e·Å^{-3} -> e·Bohr^{-3})\n', convert_Ang3_to_bohr3);

%% =========================================================================
%%  Part A: Cluster processing (spherical box, R=12 Å)
%% =========================================================================
sphere_radius = 12.0 / bohr2ang;           % au

cluster_grid_x = linspace(-sphere_radius, sphere_radius, Nx_c);
cluster_grid_y = linspace(-sphere_radius, sphere_radius, Ny_c);
cluster_grid_z = linspace(-sphere_radius, sphere_radius, Nz_c);

% Surface F⁻ atom reference position
z_offset = 0.5 * lattice_constant;
xFc = -0.8 * lattice_constant;
yFc =  0.0;
zFc =  0.0 + z_offset;

fprintf('\n=== Cluster ===\n');
fprintf('F atom theoretical coords: (%.4f, %.4f, %.4f) au\n', xFc, yFc, zFc);

[~, ixFc] = min(abs(cluster_grid_x - xFc));
[~, iyFc] = min(abs(cluster_grid_y - yFc));
[~, izFc] = min(abs(cluster_grid_z - zFc));
xFc_g = cluster_grid_x(ixFc); yFc_g = cluster_grid_y(iyFc); zFc_g = cluster_grid_z(izFc);
fprintf('Nearest grid point: (%.4f, %.4f, %.4f) au, indices [iz=%d, iy=%d, ix=%d]\n', ...
    xFc_g, yFc_g, zFc_g, izFc, iyFc, ixFc);

%% =========================================================================
%%  Part B: Slab processing (auto-detection + z-axis precision peak alignment)
%% =========================================================================
slab_grid_x = (0:Nx_s-1) * slab_spacing_au;
slab_grid_y = (0:Ny_s-1) * slab_spacing_au;
slab_grid_z = (0:Nz_s-1) * slab_spacing_au;

fprintf('\n=== Slab ===\n');
fprintf('Box size: %.3f x %.3f x %.3f au\n', ...
    slab_grid_x(end)+slab_spacing_au, slab_grid_y(end)+slab_spacing_au, slab_grid_z(end)+slab_spacing_au);

% --- z-direction atomic layer detection ---
z_density_profile = squeeze(max(max(density_slab, [], 2), [], 3));
z_density_smooth  = conv(z_density_profile, [1 1 1]/3, 'same');

peaks_z = [];
for i = 2:Nz_s-1
    if z_density_smooth(i) > z_density_smooth(i-1) && z_density_smooth(i) > z_density_smooth(i+1) ...
            && z_density_smooth(i) > 0.25 * max(z_density_smooth)
        peaks_z = [peaks_z; i];
    end
end
peaks_z = sort(peaks_z);
fprintf('Detected %d atomic layers, z indices: %s\n', length(peaks_z), mat2str(peaks_z'));

% Surface candidate layer: uppermost (max z)
surface_candidates = [min(peaks_z); max(peaks_z)];
[~, surf_idx] = min(abs(surface_candidates - 175));
z_f_s = surface_candidates(surf_idx);

% --- x-y in-plane F atom location (4-connected flood fill) ---
xy_slice = squeeze(density_slab(z_f_s, :, :));
threshold = 0.3 * max(xy_slice(:));
mask = xy_slice > threshold;

labels = zeros(Ny_s, Nx_s);
current_label = 0;
for j = 1:Ny_s
    for i = 1:Nx_s
        if mask(j,i) && labels(j,i) == 0
            current_label = current_label + 1;
            queue = [j, i];
            labels(j,i) = current_label;
            qhead = 1;
            while qhead <= size(queue,1)
                qj = queue(qhead,1); qi = queue(qhead,2); qhead = qhead + 1;
                if qj > 1 && mask(qj-1,qi) && labels(qj-1,qi)==0
                    labels(qj-1,qi)=current_label; queue=[queue; qj-1, qi]; end
                if qj < Ny_s && mask(qj+1,qi) && labels(qj+1,qi)==0
                    labels(qj+1,qi)=current_label; queue=[queue; qj+1, qi]; end
                if qi > 1 && mask(qj,qi-1) && labels(qj,qi-1)==0
                    labels(qj,qi-1)=current_label; queue=[queue; qj, qi-1]; end
                if qi < Nx_s && mask(qj,qi+1) && labels(qj,qi+1)==0
                    labels(qj,qi+1)=current_label; queue=[queue; qj, qi+1]; end
            end
        end
    end
end

cx = round(Nx_s/2);
cy = round(Ny_s/2);
center_label = labels(cy, cx);

if center_label > 0
    [yy, xx] = find(labels == center_label);
    fx = round(mean(xx));
    fy = round(mean(yy));
    fprintf('Center point (%d,%d) in connected domain %d\n', cx, cy, center_label);
else
    peaks_xy = [];
    min_pixels = 5;
    for k = 1:current_label
        [yy, xx] = find(labels == k);
        if numel(yy) < min_pixels, continue; end
        px = round(mean(xx)); py = round(mean(yy));
        peak_val = max(xy_slice(sub2ind(size(xy_slice), yy, xx)));
        peaks_xy = [peaks_xy; px, py, peak_val, k];
    end
    dist_to_center = sqrt((peaks_xy(:,1) - cx).^2 + (peaks_xy(:,2) - cy).^2);
    [~, nearest_idx] = min(dist_to_center);
    fx = peaks_xy(nearest_idx, 1);
    fy = peaks_xy(nearest_idx, 2);
    fprintf('Center in gap, selected nearest domain %d\n', peaks_xy(nearest_idx, 4));
end

% --- Critical fix: z-direction precise peak alignment ---
z_line_s = squeeze(density_slab(:, fy, fx));   % line density along z at (fx, fy)
[~, iz_peak_s] = max(z_line_s);
z_peak_s = slab_grid_z(iz_peak_s);
fprintf('Slab F layer coarse z: %.3f au, precise peak at (fx,fy) z: %.3f au\n', ...
    slab_grid_z(z_f_s), z_peak_s);

% Replace coarse with precise peak for strict alignment with Cluster zFc
z_f_s = iz_peak_s;
zFs = z_peak_s;
xFs = slab_grid_x(fx);
yFs = slab_grid_y(fy);

fprintf('Slab surface F: indices(%d,%d,%d) -> coords(%.3f, %.3f, %.3f) au\n', ...
    fx, fy, z_f_s, xFs, yFs, zFs);

%% Alignment verification
fprintf('\n=== Alignment Verification ===\n');
fprintf('Cluster F nucleus z: %.4f au\n', zFc);
fprintf('Slab   F nucleus z: %.4f au\n', zFs);
fprintf('z-direction offset: %.4f au\n', zFc - zFs);

%% =========================================================================
%%  Part C: Extract local regions on matching recorded grids
%% =========================================================================
% --- Cluster local ---
ix_c = find(cluster_grid_x >= xFc - box_half_width & cluster_grid_x <= xFc + box_half_width);
iy_c = find(cluster_grid_y >= yFc - box_half_width & cluster_grid_y <= yFc + box_half_width);
iz_c = find(cluster_grid_z >= zFc - box_half_width & cluster_grid_z <= zFc + box_half_width);

density_cluster_local = density_cluster(iz_c, iy_c, ix_c);
x_c_loc = cluster_grid_x(ix_c); y_c_loc = cluster_grid_y(iy_c); z_c_loc = cluster_grid_z(iz_c);

% --- Slab local ---
ix_s = find(slab_grid_x >= xFs - box_half_width & slab_grid_x <= xFs + box_half_width);
iy_s = find(slab_grid_y >= yFs - box_half_width & slab_grid_y <= yFs + box_half_width);
iz_s = find(slab_grid_z >= zFs - box_half_width & slab_grid_z <= zFs + box_half_width);

density_slab_local = density_slab(iz_s, iy_s, ix_s);
x_s_loc = slab_grid_x(ix_s); y_s_loc = slab_grid_y(iy_s); z_s_loc = slab_grid_z(iz_s);

% Compare only existing samples. Grids must coincide after centering.
x_uni = x_c_loc - xFc;
y_uni = y_c_loc - yFc;
z_uni = z_c_loc - zFc;
source_axes = {x_s_loc - xFs, y_s_loc - yFs, z_s_loc - zFs};
target_axes = {x_uni, y_uni, z_uni};
[density_cluster_uniform, density_slab_uniform, cell_volume] = ...
    compare_on_native_grids(density_cluster_local, density_slab_local, target_axes, source_axes);

fprintf('\nUnified grid dimensions: [%d, %d, %d]\n', length(y_uni), length(x_uni), length(z_uni));

%% =========================================================================
%%  Part D: Fixed spherical cutoff (R_cut = 2.8 au)
%% =========================================================================
fprintf('Fixed cutoff radius: %.1f au (%.3f Å)\n', cutoff_radius, cutoff_radius*bohr2ang);

[X_mg, Y_mg, Z_mg] = meshgrid(x_uni, y_uni, z_uni);   % [Ny, Nx, Nz]
radial_distance = sqrt(X_mg.^2 + Y_mg.^2 + Z_mg.^2);

% Set outside-sphere to NaN
density_cluster_uniform(radial_distance > cutoff_radius) = NaN;
density_slab_uniform(radial_distance > cutoff_radius) = NaN;

% Differential density
density_diff = density_cluster_uniform - density_slab_uniform;

%% =========================================================================
%%  Part E: Charge center coordinates (density-weighted centroid)
%% =========================================================================
fprintf('\n=== Charge Center Coordinates (density-weighted centroid) ===\n');

mask_sphere = ~isnan(density_cluster_uniform);

% --- Cluster ---
density_cluster_total = sum(density_cluster_uniform(mask_sphere), 'omitnan');
x_cm_c = sum(density_cluster_uniform(mask_sphere) .* X_mg(mask_sphere), 'omitnan') / density_cluster_total;
y_cm_c = sum(density_cluster_uniform(mask_sphere) .* Y_mg(mask_sphere), 'omitnan') / density_cluster_total;
z_cm_c = sum(density_cluster_uniform(mask_sphere) .* Z_mg(mask_sphere), 'omitnan') / density_cluster_total;

fprintf('Cluster F^- charge center:\n');
fprintf('  Relative coords: (%.5f, %.5f, %.5f) au\n', x_cm_c, y_cm_c, z_cm_c);

% --- Slab ---
density_slab_total = sum(density_slab_uniform(mask_sphere), 'omitnan');
x_cm_s = sum(density_slab_uniform(mask_sphere) .* X_mg(mask_sphere), 'omitnan') / density_slab_total;
y_cm_s = sum(density_slab_uniform(mask_sphere) .* Y_mg(mask_sphere), 'omitnan') / density_slab_total;
z_cm_s = sum(density_slab_uniform(mask_sphere) .* Z_mg(mask_sphere), 'omitnan') / density_slab_total;

fprintf('Slab F^- charge center:\n');
fprintf('  Relative coords: (%.5f, %.5f, %.5f) au\n', x_cm_s, y_cm_s, z_cm_s);

dx = x_cm_c - x_cm_s; dy = y_cm_c - y_cm_s; dz = z_cm_c - z_cm_s;
fprintf('Charge center shift (Cluster - Slab): (%.5f, %.5f, %.5f) au\n', dx, dy, dz);
fprintf('Shift magnitude: %.5f au (%.4f Å)\n', sqrt(dx^2+dy^2+dz^2), sqrt(dx^2+dy^2+dz^2)*bohr2ang);

%% =========================================================================
%%  Part F: 2D slice (z = 0) — Figures (a)(b)
%% =========================================================================
[~, iz0] = min(abs(z_uni));

density_cluster_xy = squeeze(density_cluster_uniform(:, :, iz0));
density_slab_xy    = squeeze(density_slab_uniform(:, :, iz0));

R_xy = sqrt(X_mg(:,:,iz0).^2 + Y_mg(:,:,iz0).^2);
density_cluster_xy(R_xy > cutoff_radius) = NaN;
density_slab_xy(R_xy > cutoff_radius)    = NaN;

%% =========================================================================
%%  Part G: Radial analysis — Figure C: Δρ(r)
%% =========================================================================
radial_step = 0.1;
r_edges = 0 : radial_step : cutoff_radius;
r_centers = r_edges(1:end-1) + radial_step/2;

density_cluster_radial = zeros(size(r_centers));
density_slab_radial    = zeros(size(r_centers));

for i_r = 1:length(r_centers)
    mask_radial = radial_distance >= r_edges(i_r) & radial_distance < r_edges(i_r+1);
    if nnz(mask_radial) > 0
        density_cluster_radial(i_r) = mean(density_cluster_uniform(mask_radial), 'omitnan');
        density_slab_radial(i_r)    = mean(density_slab_uniform(mask_radial), 'omitnan');
    else
        density_cluster_radial(i_r) = NaN;
        density_slab_radial(i_r)    = NaN;
    end
end

density_diff_radial = density_cluster_radial - density_slab_radial;

fprintf('\n=== Radial Analysis ===\n');
[~, max_idx] = max(density_diff_radial);
[~, min_idx] = min(density_diff_radial);
fprintf('Δρ(r) max positive peak: r = %.2f au, Δρ = %+.4e\n', r_centers(max_idx), density_diff_radial(max_idx));
fprintf('Δρ(r) max negative peak: r = %.2f au, Δρ = %+.4e\n', r_centers(min_idx), density_diff_radial(min_idx));

%% =========================================================================
%%  Part H: z-axis density difference — Figure D: Δρ(z) (x=y=0)
%% =========================================================================
[~, ix0] = min(abs(x_uni));   % x = 0
[~, iy0] = min(abs(y_uni));   % y = 0

density_diff_z = squeeze(density_diff(iy0, ix0, :));   % along z-axis

% Diagnostics: +z vs -z hemisphere integrals
mask_z_pos = z_uni > 0 & z_uni <= cutoff_radius;
mask_z_neg = z_uni < 0 & abs(z_uni) <= cutoff_radius;

if nnz(mask_z_pos) > 0 && nnz(mask_z_neg) > 0
    delta_Q_pos = sum(density_diff_z(mask_z_pos), 'omitnan') * cell_volume;
    delta_Q_neg = sum(density_diff_z(mask_z_neg), 'omitnan') * cell_volume;
    fprintf('\n=== z-axis Difference Diagnostics ===\n');
    fprintf('+z hemisphere (vacuum side) net charge diff: %+.4e e\n', delta_Q_pos);
    fprintf('-z hemisphere (substrate side) net charge diff: %+.4e e\n', delta_Q_neg);
    if delta_Q_pos > abs(delta_Q_neg)
        fprintf('-> +z direction net charge surplus, confirming electron density spill-out toward vacuum\n');
    elseif delta_Q_pos > 0
        fprintf('-> +z direction has net charge surplus, supporting vacuum-side delocalization\n');
    end
end

%% =========================================================================
%%  Part I: Write output data files
%% =========================================================================
fprintf('\n=== Saving Data Files ===\n');

% (a) Cluster 2D
fid = fopen('cluster_xy_z0.txt', 'w');
fprintf(fid, 'x(au)\ty(au)\trho(e·Bohr^{-3})\n');
for j = 1:length(y_uni)
    for i = 1:length(x_uni)
        if isnan(density_cluster_xy(j,i))
            fprintf(fid, '%.6f\t%.6f\tNaN\n', x_uni(i), y_uni(j));
        else
            fprintf(fid, '%.6f\t%.6f\t%.6e\n', x_uni(i), y_uni(j), density_cluster_xy(j,i));
        end
    end
end
fclose(fid);
fprintf('Saved: cluster_xy_z0.txt\n');

% (b) Slab 2D
fid = fopen('slab_xy_z0.txt', 'w');
fprintf(fid, 'x(au)\ty(au)\trho(e·Bohr^{-3})\n');
for j = 1:length(y_uni)
    for i = 1:length(x_uni)
        if isnan(density_slab_xy(j,i))
            fprintf(fid, '%.6f\t%.6f\tNaN\n', x_uni(i), y_uni(j));
        else
            fprintf(fid, '%.6f\t%.6f\t%.6e\n', x_uni(i), y_uni(j), density_slab_xy(j,i));
        end
    end
end
fclose(fid);
fprintf('Saved: slab_xy_z0.txt\n');

% (c) Figure C: radial density difference Δρ(r)
fid = fopen('radial_density_diff.txt', 'w');
fprintf(fid, 'r(au)\tΔρ(e·Bohr^{-3})\n');
for i = 1:length(r_centers)
    fprintf(fid, '%.6f\t%.6e\n', r_centers(i), density_diff_radial(i));
end
fclose(fid);
fprintf('Saved: radial_density_diff.txt (Figure C: Δρ(r))\n');

% (d) Figure D: z-axis density difference Δρ(z)
fid = fopen('z_axis_density_diff.txt', 'w');
fprintf(fid, 'z(au)\tΔρ(e·Bohr^{-3})\n');
for i = 1:length(z_uni)
    if ~isnan(density_diff_z(i))
        fprintf(fid, '%.6f\t%.6e\n', z_uni(i), density_diff_z(i));
    end
end
fclose(fid);
fprintf('Saved: z_axis_density_diff.txt (Figure D: Δρ(z))\n');

% (e) Charge center summary
fid = fopen('charge_center_summary.txt', 'w');
fprintf(fid, '=== LiF Surface F^- Analysis ===\n');
fprintf(fid, 'Cutoff: %.1f au (%.3f Å)\n', cutoff_radius, cutoff_radius*bohr2ang);
fprintf(fid, 'Cluster F nucleus: (%.6f, %.6f, %.6f) au\n', xFc, yFc, zFc);
fprintf(fid, 'Slab   F nucleus: (%.6f, %.6f, %.6f) au\n', xFs, yFs, zFs);
fprintf(fid, 'z alignment offset: %.6f au\n', zFc - zFs);
fprintf(fid, '\nCluster charge center: (%.6f, %.6f, %.6f) au\n', x_cm_c, y_cm_c, z_cm_c);
fprintf(fid, 'Slab   charge center: (%.6f, %.6f, %.6f) au\n', x_cm_s, y_cm_s, z_cm_s);
fprintf(fid, 'Shift (C-S): (%.6f, %.6f, %.6f) au\n', dx, dy, dz);
fprintf(fid, 'Shift magnitude: %.6f au = %.6f Å\n', ...
    sqrt(dx^2+dy^2+dz^2), sqrt(dx^2+dy^2+dz^2)*bohr2ang);
fclose(fid);
fprintf('Saved: charge_center_summary.txt\n');

%% ---- z-axis cumulative difference ----
density_diff_z_step = z_uni(2) - z_uni(1);

% +z hemisphere (vacuum side) cumulative
delta_Q_z_pos = sum(density_diff_z(mask_z_pos), 'omitnan') * cell_volume;  % e

% -z hemisphere (substrate side) cumulative
delta_Q_z_neg = sum(density_diff_z(mask_z_neg), 'omitnan') * cell_volume;  % e

fprintf('+z hemisphere (vacuum side) net charge diff: %+.4e e\n', delta_Q_z_pos);
fprintf('-z hemisphere (substrate side) net charge diff: %+.4e e\n', delta_Q_z_neg);

% Cumulative curve
delta_Q_z_cum = cumsum(density_diff_z(~isnan(density_diff_z))) * cell_volume;
z_valid = z_uni(~isnan(density_diff_z));

% Output cumulative data
fid = fopen('z_axis_cumulative.txt', 'w');
fprintf(fid, 'z(au)\tΔQ_z(e)\n');
for i = 1:length(z_valid)
    fprintf(fid, '%.6f\t%.6e\n', z_valid(i), delta_Q_z_cum(i));
end
fclose(fid);

fprintf('\n=== All Complete ===\n');
fprintf('File list:\n');
fprintf('  1. cluster_xy_z0.txt         - (a) Cluster z=0 slice\n');
fprintf('  2. slab_xy_z0.txt            - (b) Slab z=0 slice\n');
fprintf('  3. radial_density_diff.txt   - (c) Δρ(r) [e·Bohr^{-3}]\n');
fprintf('  4. z_axis_density_diff.txt   - (d) Δρ(z) [e·Bohr^{-3}]\n');
fprintf('  5. charge_center_summary.txt - Charge center & alignment info\n');
