%% =============================================================================
%  calc_polarizability - Calculate polarizability tensor components from
%                        finite-field DFT calculations (Octopus)
%
%  This script reads 3D electron density from three Octopus static DFT
%  calculations (zero field, +E_z, -E_z) and computes the polarizability
%  tensor components α_xz, α_yz, α_zz for a surface F⁻ ion:
%
%    α_ij = (μ_i(+E_z) - μ_i(-E_z)) / (2 · E_z)
%
%  where μ_i = -∫ ρ(r) · r_i d³r (electronic dipole moment).
%
%  The script automatically:
%    1. Parses Octopus inp file for grid parameters and atom coordinates
%    2. Locates the target F⁻ ion (auto/manual/density-minimum modes)
%    3. Reads NetCDF density files with correct dimension ordering (z,y,x)
%    4. Handles unit conversion (eV_Angstrom → atomic units)
%    5. Constructs offset coordinate grid centered on F⁻ ion
%    6. Computes dipole moments and polarizability within spherical cutoff
%
%  OUTPUT (console):
%    - Zero-field dipole moments μ_x(0), μ_y(0), μ_z(0) [e·au]
%    - Finite-field dipole moments μ_i(±E_z)
%    - Polarizability components α_xz, α_yz, α_zz [au³] and [Å³]
%    - Linearity check
%
%  USAGE:
%    1. Set CONFIGURATION paths for zero-field, +E, -E density files
%    2. Set inp_dir to Octopus input file directory
%    3. Adjust R_cut, E_field, and center_mode as needed
%    4. Octave users: ensure netcdf package is loaded:
%         pkg install -forge netcdf
%         pkg load netcdf
%       Or run: run('matlab/utils/octave_startup.m');
%    5. Run the script
%
%  DEPENDENCIES:
%    - MATLAB / GNU Octave 4.4+
%    - Octopus inp file (for grid parameters and atom positions)
%    - Octopus NetCDF density files (static/density.ncdf)
%    - parse_octopus_inp_v7(), eval_expr_with_vars_v2(), parse_coord_lines_v7(),
%      eval_coord_expr_v7(), parse_value_with_unit(), parse_vector_with_unit(),
%      convert_to_au() — all defined as subfunctions in this file
%
%  Octave USERS: ensure netcdf package is loaded (see octave_startup.m)
%
%  REFERENCES:
%    Project: "Occupied-Space Constraints and Finite-Time Recovery in
%             Electron Detachment during F-/LiF(100) Scattering"
%    Code: Octopus DFT (https://octopus-code.org)
%    Method: Finite-field polarizability, see Umari & Pasquarello, PRL 89, 157602 (2002)
% =============================================================================
clear; clc; close all;

%% ==================== CONFIGURATION ====================
% Density file paths for zero-field, +E, and -E calculations (USER: adjust)
cube_zero = "Z:\version1.5\LiF\v1.3\\static\\density.ncdf";
cube_plus = "Z:\version1.5\LiF\v1.3\\polarize+\\static\\density.ncdf";
cube_minus = "Z:\version1.5\LiF\v1.3\\polarize-\\static\\density.ncdf";

% Octopus input directory (USER: adjust)
inp_dir = "Z:\version1.5\LiF\v1.3\\";

% Analysis parameters
manual_spacing  = 0.1890;   % au, fallback grid spacing
E_field         = 0.002;    % au, applied electric field (z-direction)
cutoff_radius   = 3.0;      % au, spherical integration cutoff

% F⁻ ion location mode (USER: choose 'auto', 'manual', or 'density')
center_mode = 'auto';
manual_center_idx = [131, 121, 105];   % [iz, iy, ix] when center_mode='manual'
f_selection = 'center';      % 'center', 'highest_z', or 'lowest_z'

% Manual grid origin override (leave empty for auto-detection)
origin_manual = [];          % [x0, y0, z0] au
% ========================================================

%% ---- Locate Octopus input file ----
fprintf("=== Locating Octopus inp file ===\n");
if isempty(inp_dir), inp_dir = fileparts(cube_zero); end

inp_candidates = {fullfile(inp_dir, "inp"), ...
                  fullfile(inp_dir, "inp.inp"), ...
                  fullfile(inp_dir, "octopus.inp"), ...
                  fullfile(inp_dir, "input"), ...
                  fullfile(inp_dir, "INPUT")};

inp_path = "";
for i = 1:length(inp_candidates)
    if isfile(inp_candidates{i})
        inp_path = inp_candidates{i};
        fprintf("  Found: %s\n", inp_path);
        break;
    end
end
if inp_path == ""
    fprintf("  Warning: No standard inp file found, using manual fallback parameters...\n");
    use_manual_inp = true;
else
    use_manual_inp = false;
end

%% ---- Parse inp parameters & variables ----
fprintf("=== Reading Octopus input parameters ===\n");
if ~use_manual_inp
    [vars, inp] = parse_octopus_inp_v7(inp_path);
else
    vars = struct();
    inp = struct('UnitsInput','atomic','UnitsOutput','atomic',...
                 'Spacing',[],'SpacingUnit','','BoxShape','parallelepiped',...
                 'Lsize',[],'Radius',[],'Atoms',struct('Symbol',{},'Coords',{}));
end
fprintf("  UnitsInput  : %s\n", inp.UnitsInput);
fprintf("  UnitsOutput : %s\n", inp.UnitsOutput);
fprintf("  BoxShape    : %s\n", inp.BoxShape);

%% ---- Read NetCDF density, explicit dimension order (z, y, x) ----
fprintf("\n=== Reading NetCDF density files ===\n");
density_zero  = ncread(cube_zero,  "rdata");
density_plus  = ncread(cube_plus,  "rdata");
density_minus = ncread(cube_minus, "rdata");

% Dimension order: dim 1 = z, dim 2 = y, dim 3 = x
nz = size(density_zero, 1);
ny = size(density_zero, 2);
nx = size(density_zero, 3);
fprintf("  Grid dimensions (z,y,x): %d x %d x %d\n", nz, ny, nx);

% Try to read NetCDF coordinate variables (for validation / priority use)
try
    x_nc = ncread(cube_zero, 'x');
    y_nc = ncread(cube_zero, 'y');
    z_nc = ncread(cube_zero, 'z');
    nc_has_coords = true;
catch
    nc_has_coords = false;
    fprintf("  No coordinate variables x/y/z in NetCDF, relying on inp or manual fallback.\n");
end

if nc_has_coords
    if length(x_nc) ~= nx || length(y_nc) ~= ny || length(z_nc) ~= nz
        warning('Coordinate variable lengths do not match density dimensions! Check file.');
    end
    nc_x = x_nc(:)';
    nc_y = y_nc(:)';
    nc_z = z_nc(:)';
    fprintf("  NetCDF coordinate ranges:\n");
    fprintf("    x: [%.4f, %.4f], y: [%.4f, %.4f], z: [%.4f, %.4f]\n", ...
        nc_x(1), nc_x(end), nc_y(1), nc_y(end), nc_z(1), nc_z(end));
else
    nc_x = []; nc_y = []; nc_z = [];
end

%% ---- Output unit system determination ----
fprintf("\n=== Output unit system determination ===\n");
units_output = lower(inp.UnitsOutput);
is_angstrom_output = contains(units_output, 'angstrom') || contains(units_output, 'ev_angstrom');
if is_angstrom_output
    fprintf("  UnitsOutput = %s: Octopus density output in e/Angstrom^3\n", inp.UnitsOutput);
    fprintf("  Note: inp internal parameters still in atomic units, only output units changed\n");
    if nc_has_coords
        nc_x = nc_x / 0.529177210903;
        nc_y = nc_y / 0.529177210903;
        nc_z = nc_z / 0.529177210903;
        fprintf("  NetCDF coordinates converted from Angstrom to au\n");
    end
else
    fprintf("  UnitsOutput = atomic: density output in e/au^3\n");
end

%% ---- Determine grid spacing (dx, dy, dz) in au ----
fprintf("\n=== Determining grid spacing ===\n");
if ~isempty(inp.Spacing) && isfinite(inp.Spacing(1))
    dx_inp = inp.Spacing(1);
    if ~isempty(inp.SpacingUnit)
        dx_inp = convert_to_au(dx_inp, inp.SpacingUnit);
        fprintf("  Spacing (inp, explicit unit %s): %.4f au\n", inp.SpacingUnit, dx_inp);
    else
        if contains(inp.UnitsInput, 'angstrom') || contains(inp.UnitsInput, 'ev')
            dx_inp = dx_inp * 1.8897261339;
            fprintf("  Spacing (inp, UnitsInput=%s): %.4f au\n", inp.UnitsInput, dx_inp);
        else
            fprintf("  Spacing (inp, atomic): %.4f au\n", dx_inp);
        end
    end
    dx = dx_inp; dy = dx_inp; dz = dx_inp;
    if length(inp.Spacing) >= 3
        dy = convert_to_au(inp.Spacing(2), inp.SpacingUnit);
        dz = convert_to_au(inp.Spacing(3), inp.SpacingUnit);
    end
elseif nc_has_coords && length(nc_x) > 1
    dx = abs(nc_x(2) - nc_x(1));
    dy = abs(nc_y(2) - nc_y(1));
    dz = abs(nc_z(2) - nc_z(1));
    fprintf("  Spacing inferred from NetCDF coords: dx=%.4f, dy=%.4f, dz=%.4f au\n", dx, dy, dz);
else
    dx = manual_spacing; dy = manual_spacing; dz = manual_spacing;
    fprintf("  Using manual fallback Spacing: %.4f au\n", dx);
end
dV = dx * dy * dz;
fprintf("  Final dV = %.6f au^3\n", dV);

%% ---- Determine grid origin (x0, y0, z0) in au ----
fprintf("\n=== Determining grid origin ===\n");
if ~isempty(origin_manual) && length(origin_manual) == 3
    origin = origin_manual(:)';  % [x0, y0, z0]
    fprintf("  Using manual origin: [%.4f, %.4f, %.4f] au\n", origin);
elseif nc_has_coords && ~isempty(nc_x)
    origin = [nc_x(1), nc_y(1), nc_z(1)];
    fprintf("  Origin from NetCDF coords: [%.4f, %.4f, %.4f] au\n", origin);
elseif contains(inp.BoxShape, "parallelepiped") || isempty(inp.BoxShape)
    if ~isempty(inp.Lsize) && length(inp.Lsize) >= 3
        L = inp.Lsize(1:3);
    elseif ~isempty(inp.Lsize) && isfinite(inp.Lsize(1))
        L = [inp.Lsize(1), inp.Lsize(1), inp.Lsize(1)];
    else
        L = [0,0,0];
    end
    origin = -L/2;
    fprintf("  parallelepiped: origin = -Lsize/2 = [%.4f, %.4f, %.4f] au\n", origin);
elseif contains(inp.BoxShape, "sphere")
    if ~isempty(inp.Radius) && isfinite(inp.Radius)
        R_box = inp.Radius;
        origin = [-R_box, -R_box, -R_box];
        fprintf("  sphere mode: origin = -Radius = [%.4f, %.4f, %.4f] au\n", origin);
        expected_range = 2 * R_box;
        actual_range_x = (nx-1) * dx;
        if abs(expected_range - actual_range_x) > 0.01*dx
            fprintf("  Warning: Radius inconsistent with grid x range!\n");
        end
    else
        origin = -[(nx-1)/2*dx, (ny-1)/2*dy, (nz-1)/2*dz];
        fprintf("  sphere mode: auto-inferred origin = [%.4f, %.4f, %.4f] au\n", origin);
    end
else
    origin = [0, 0, 0];
    warning("Unknown BoxShape, origin set to [0,0,0]");
end

%% ---- Locate target F⁻ ion ----
fprintf("\n=== Locating F⁻ ion ===\n");
switch lower(center_mode)
    case 'auto'
        if isempty(inp.Atoms) || isempty(inp.Atoms.Coords)
            error("inp contains no atom coordinates! Switch to center_mode='density' or 'manual'");
        end
        is_F = strcmpi(inp.Atoms.Symbol, 'F');
        if ~any(is_F)
            fprintf("  No pure F found, trying symbols containing 'F' but not 'udf'...\n");
            is_F = contains(lower(inp.Atoms.Symbol), 'f') & ~contains(lower(inp.Atoms.Symbol), 'udf');
            if ~any(is_F)
                error("No F atoms found!");
            end
        end
        F_coords = inp.Atoms.Coords(is_F, :);  % each row (x, y, z)
        fprintf("  Found %d pure F atoms:\n", size(F_coords,1));
        for i = 1:size(F_coords,1)
            fprintf("    F(%d): (%.4f, %.4f, %.4f) au\n", i, F_coords(i,:));
        end

        % Boundary check
        grid_min = origin;
        grid_max = origin + [(nx-1)*dx, (ny-1)*dy, (nz-1)*dz];
        for i = 1:size(F_coords,1)
            for dim = 1:3
                if F_coords(i,dim) < grid_min(dim) || F_coords(i,dim) > grid_max(dim)
                    fprintf("  Warning: F(%d) coord %.4f au outside grid range!\n", i, F_coords(i,dim));
                end
            end
        end

        if size(F_coords,1) == 1
            target_F = F_coords(1,:);
        else
            switch lower(f_selection)
                case 'center'
                    gc = mean(inp.Atoms.Coords,1);
                    [~,idx] = min(vecnorm(F_coords - gc, 2, 2));
                    fprintf("  Selected F(%d) nearest geometric center\n", idx);
                case 'highest_z'
                    [~,idx] = max(F_coords(:,3));
                case 'lowest_z'
                    [~,idx] = min(F_coords(:,3));
                otherwise
                    idx = 1;
            end
            target_F = F_coords(idx,:);
        end

        % Physical coords to indices (idx_x, idx_y, idx_z) -> stored as [iz, iy, ix]
        idx_x = round((target_F(1) - origin(1)) / dx) + 1;
        idx_y = round((target_F(2) - origin(2)) / dy) + 1;
        idx_z = round((target_F(3) - origin(3)) / dz) + 1;
        idx_x = max(min(idx_x, nx), 1);
        idx_y = max(min(idx_y, ny), 1);
        idx_z = max(min(idx_z, nz), 1);
        idx_center = [idx_z, idx_y, idx_x];   % matches rho dimension order (z,y,x)
        fprintf("  Target coords: (%.4f, %.4f, %.4f) au\n", target_F);
        fprintf("  Grid indices (z,y,x): [%d, %d, %d]\n", idx_center);

    case 'manual'
        idx_center = manual_center_idx;  % user-provided [iz, iy, ix]
        fprintf("  Manual indices (z,y,x): [%d, %d, %d]\n", idx_center);

    case 'density'
        border = round(0.05 * [nz, ny, nx]);
        density_interior = density_zero(border(1)+1:nz-border(1), border(2)+1:ny-border(2), border(3)+1:nx-border(3));
        [min_val, min_idx] = min(density_interior(:));
        [izr, iyr, ixr] = ind2sub(size(density_interior), min_idx);
        idx_center = [izr+border(1), iyr+border(2), ixr+border(3)];
        fprintf("  Density minimum %.4f at indices (z,y,x): [%d, %d, %d]\n", min_val, idx_center);

    otherwise
        error("Unknown center_mode: %s", center_mode);
end

%% ---- Build coordinate grid & spherical mask ----
% Offset grid (centered on target ion)
z0 = ((1:nz) - idx_center(1)) * dz;
y0 = ((1:ny) - idx_center(2)) * dy;
x0 = ((1:nx) - idx_center(3)) * dx;
[Z, Y, X] = ndgrid(z0, y0, x0);   % dim1=Z, dim2=Y, dim3=X
R = sqrt(X.^2 + Y.^2 + Z.^2);
mask = R < cutoff_radius;

%% ---- Density unit conversion (eV_Angstrom -> au) ----
if is_angstrom_output
    au3_per_A3 = 0.148184;   % 0.529177^3
    density_zero  = density_zero  * au3_per_A3;
    density_plus  = density_plus  * au3_per_A3;
    density_minus = density_minus * au3_per_A3;
    fprintf("\n  [eV_Angstrom] Density converted from e/A^3 to e/au^3\n");
end

%% ---- Diagnostics and dipole calculation (all three directions) ----
electron_count_zero = sum(density_zero(mask), 'all') * dV;
fprintf("\nElectrons within R_cut=%.1f au: %.4f\n", cutoff_radius, electron_count_zero);

N_total = sum(density_zero(:), 'all') * dV;
fprintf("Full grid integrated electron count: %.4f e\n", N_total);
fprintf("Density range: [%.4f, %.4f] e/au^3\n", min(density_zero(:)), max(density_zero(:)));

if electron_count_zero < 1e-6
    error("Electron count in sphere near 0! Check F⁻ index and R_cut");
end

% Three-direction dipole moments
mu_x_p = -sum(density_plus(mask)  .* X(mask), 'all') * dV;
mu_x_m = -sum(density_minus(mask) .* X(mask), 'all') * dV;
mu_x_0 = -sum(density_zero(mask)  .* X(mask), 'all') * dV;

mu_y_p = -sum(density_plus(mask)  .* Y(mask), 'all') * dV;
mu_y_m = -sum(density_minus(mask) .* Y(mask), 'all') * dV;
mu_y_0 = -sum(density_zero(mask)  .* Y(mask), 'all') * dV;

mu_z_p = -sum(density_plus(mask)  .* Z(mask), 'all') * dV;
mu_z_m = -sum(density_minus(mask) .* Z(mask), 'all') * dV;
mu_z_0 = -sum(density_zero(mask)  .* Z(mask), 'all') * dV;

% Polarizability components (E_z ≠ 0 only)
alpha_xz = (mu_x_p - mu_x_m) / (2.0 * E_field);
alpha_yz = (mu_y_p - mu_y_m) / (2.0 * E_field);
alpha_zz = (mu_z_p - mu_z_m) / (2.0 * E_field);

%% ---- Results output ----
fprintf("\n=== Zero-field dipole moments (e·au) ===\n");
fprintf("μ_x(0) = %.6f\n", mu_x_0);
fprintf("μ_y(0) = %.6f\n", mu_y_0);
fprintf("μ_z(0) = %.6f\n", mu_z_0);

fprintf("\n=== Finite-field dipole moments (E_z=%.4f a.u.) ===\n", E_field);
fprintf("μ_x(+E) = %.6f, μ_x(-E) = %.6f\n", mu_x_p, mu_x_m);
fprintf("μ_y(+E) = %.6f, μ_y(-E) = %.6f\n", mu_y_p, mu_y_m);
fprintf("μ_z(+E) = %.6f, μ_z(-E) = %.6f\n", mu_z_p, mu_z_m);

fprintf("\n=== Polarizability components (au³) [E_z direction field only] ===\n");
fprintf("α_xz = %.4f au³ = %.4f Å³\n", alpha_xz, alpha_xz * 0.148184);
fprintf("α_yz = %.4f au³ = %.4f Å³\n", alpha_yz, alpha_yz * 0.148184);
fprintf("α_zz = %.4f au³ = %.4f Å³\n", alpha_zz, alpha_zz * 0.148184);

% Linearity check (based on z component)
mu_ind_p_z = mu_z_p - mu_z_0;
alpha_zz_alt = mu_ind_p_z / E_field;
if abs(alpha_zz) > 1e-10
    rel_diff = abs(alpha_zz - alpha_zz_alt) / abs(alpha_zz);
    fprintf("\nα_zz linearity difference: %.2f%%\n", rel_diff*100);
    if rel_diff < 0.05
        fprintf("Linear response confirmed.\n");
    else
        fprintf("Warning: Significant nonlinearity, consider reducing E_field.\n");
    end
end

fprintf("\nNote: α_xz, α_yz are cross-terms; non-zero values indicate non-centrosymmetric F⁻ environment.\n");
fprintf("To obtain α_xx, α_yy, apply E_x and E_y fields respectively.\n");

%% ########################################################################
%% ########################### SUBFUNCTIONS ################################
%% ########################################################################

function [vars, inp] = parse_octopus_inp_v7(filepath)
    vars = struct();
    inp = struct('UnitsInput','atomic','UnitsOutput','atomic',...
                 'Spacing',[],'SpacingUnit','','BoxShape','',...
                 'Lsize',[],'Radius',[],'Atoms',struct('Symbol',{},'Coords',{}));
    fid = fopen(filepath, 'r');
    if fid == -1, warning("Cannot open: %s", filepath); return; end
    in_coords = false;
    coord_lines = {};
    known_params = {'spacing','unitsinput','units','unitsoutput','boxshape','lsize','radius'};
    while ~feof(fid)
        line = fgetl(fid);
        if ~ischar(line), break; end
        line_clean = regexprep(line, '[#!].*$', '');
        line_clean = strtrim(line_clean);
        if isempty(line_clean), continue; end
        if ~in_coords
            if length(line_clean) >= 12 && strcmpi(line_clean(1:12), '%Coordinates')
                in_coords = true; continue;
            elseif strcmpi(line_clean, 'Coordinates')
                in_coords = true; continue;
            end
            tokens = regexp(line_clean, '^(\w+)\s*=\s*(.+)$', 'tokens');
            if ~isempty(tokens)
                var_name = lower(strtrim(tokens{1}{1}));
                var_val = strtrim(tokens{1}{2});
                if ismember(var_name, known_params)
                    switch var_name
                        case 'spacing'
                            [num, unit] = parse_value_with_unit(var_val);
                            inp.Spacing = num; inp.SpacingUnit = unit;
                        case {'unitsinput','units'}
                            inp.UnitsInput = lower(var_val);
                        case 'unitsoutput'
                            inp.UnitsOutput = lower(var_val);
                        case 'boxshape'
                            inp.BoxShape = lower(var_val);
                        case 'lsize'
                            inp.Lsize = parse_vector_with_unit(var_val);
                        case 'radius'
                            [num, unit] = parse_value_with_unit(var_val);
                            inp.Radius = convert_to_au(num, unit);
                    end
                else
                    num_val = eval_expr_with_vars_v2(var_val, vars);
                    if ~isnan(num_val), vars.(var_name) = num_val; end
                end
            end
        else
            tmp = strtrim(line_clean);
            if strcmpi(tmp, '%') || strcmpi(tmp, 'end') || ...
                    (length(tmp) >= 1 && tmp(1) == '%')
                in_coords = false; continue;
            end
            coord_lines{end+1} = line_clean;
        end
    end
    fclose(fid);
    if ~isempty(coord_lines)
        inp.Atoms = parse_coord_lines_v7(coord_lines, inp.UnitsInput, vars);
    end
end

function val = eval_expr_with_vars_v2(expr, vars)
    expr = strtrim(expr);
    if isempty(expr), val = NaN; return; end
    if ~isempty(vars) && isstruct(vars)
        var_names = fieldnames(vars);
        [~, sort_idx] = sort(cellfun(@length, var_names), 'descend');
        var_names = var_names(sort_idx);
        for i = 1:length(var_names)
            vname = var_names{i};
            vval = vars.(vname);
            rep = num2str(vval, '%.15g');
            expr = strrep(expr, vname, rep);
            expr = strrep(expr, upper(vname), rep);
            if length(vname) > 1
                expr = strrep(expr, [upper(vname(1)), vname(2:end)], rep);
            end
        end
    end
    expr = strrep(expr, '*angstrom', '*1.8897261339');
    expr = strrep(expr, '*Angstrom', '*1.8897261339');
    expr = strrep(expr, '*ANGSTROM', '*1.8897261339');
    expr = strrep(expr, '*bohr', '*1');
    expr = strrep(expr, '*Bohr', '*1');
    expr = strrep(expr, '*au', '*1');
    expr = strrep(expr, '*AU', '*1');
    try
        val = str2num(expr);
        if isempty(val), val = NaN; end
    catch
        val = NaN;
    end
end

function atoms = parse_coord_lines_v7(lines, units_input, vars)
    symbols = {};
    coords = [];
    for i = 1:length(lines)
        line = strtrim(lines{i});
        if isempty(line), continue; end
        parts = strsplit(line, '|');
        if length(parts) < 4, continue; end
        sym = strtrim(parts{1});
        sym = regexprep(sym, '^[''"]+', '');
        sym = regexprep(sym, '[''"]+$', '');
        sym = strtrim(sym);
        if isempty(sym), continue; end
        x = eval_coord_expr_v7(strtrim(parts{2}), vars, units_input);
        y = eval_coord_expr_v7(strtrim(parts{3}), vars, units_input);
        z = eval_coord_expr_v7(strtrim(parts{4}), vars, units_input);
        if any(isnan([x,y,z]))
            fprintf("  Skipping invalid coordinate line: %s\n", line);
            continue;
        end
        symbols{end+1} = sym;
        coords = [coords; x, y, z];
    end
    atoms = struct('Symbol', {symbols}, 'Coords', coords);
end

function val = eval_coord_expr_v7(expr, vars, units_input)
    expr = strtrim(expr);
    if isempty(expr), val = NaN; return; end
    if ~isempty(vars) && isstruct(vars)
        var_names = fieldnames(vars);
        [~, sort_idx] = sort(cellfun(@length, var_names), 'descend');
        var_names = var_names(sort_idx);
        for i = 1:length(var_names)
            vname = var_names{i};
            vval = vars.(vname);
            rep = num2str(vval, '%.15g');
            expr = strrep(expr, vname, rep);
            expr = strrep(expr, upper(vname), rep);
            if length(vname) > 1
                expr = strrep(expr, [upper(vname(1)), vname(2:end)], rep);
            end
        end
    end
    expr = strrep(expr, '*angstrom', '*1.8897261339');
    expr = strrep(expr, '*Angstrom', '*1.8897261339');
    expr = strrep(expr, '*ANGSTROM', '*1.8897261339');
    expr = strrep(expr, '*bohr', '*1');
    expr = strrep(expr, '*Bohr', '*1');
    expr = strrep(expr, '*au', '*1');
    expr = strrep(expr, '*AU', '*1');
    try
        val = str2num(expr);
        if isempty(val), val = NaN; end
    catch
        val = NaN;
    end
    if contains(units_input, 'angstrom') || contains(units_input, 'ev')
        val = val * 1.8897261339;
    end
end

function [num, unit] = parse_value_with_unit(str)
    str = strtrim(str);
    tokens = regexp(str, '^(.*?)\s*\*\s*(\w+)$', 'tokens');
    if ~isempty(tokens)
        num_str = strtrim(tokens{1}{1});
        unit = lower(strtrim(tokens{1}{2}));
        num = sscanf(num_str, '%f');
        if isempty(num), num = []; end
    else
        num = sscanf(str, '%f');
        unit = '';
    end
end

function vec = parse_vector_with_unit(str)
    parts = strsplit(strtrim(str), '|');
    vec = zeros(1, length(parts));
    for i = 1:length(parts)
        [n, u] = parse_value_with_unit(strtrim(parts{i}));
        vec(i) = convert_to_au(n, u);
    end
end

function au_val = convert_to_au(val, unit)
    if isempty(unit) || isempty(val), au_val = val; return; end
    unit = lower(strtrim(unit));
    switch unit
        case {'angstrom', 'ang', 'a'}, au_val = val * 1.889726133921252;
        case {'bohr', 'au', 'a.u.', 'atomic'}, au_val = val;
        case {'nm'}, au_val = val * 18.89726133921252;
        case {'pm'}, au_val = val * 0.01889726133921252;
        otherwise, warning("Unknown unit: %s, assuming au", unit); au_val = val;
    end
end
