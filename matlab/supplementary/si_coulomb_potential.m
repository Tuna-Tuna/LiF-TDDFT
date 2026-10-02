% =========================================================================
% SI_COULOMB_POTENTIAL  Compute Coulomb potential energy map above a LiF
%                       ionic crystal surface.
%
%   This script calculates the electrostatic potential energy experienced
%   by a point-charge probe (e.g., an electron) at positions (x, y) and
%   heights z above a discrete LiF crystal slab. It serves as Supporting
%   Information for the LiF scattering study:
%
%     "Occupied-Space Constraints and Finite-Time Recovery in
%      Electron Detachment during F-/LiF(100) Scattering"
%
%   Physics context:
%     The LiF crystal is modelled as a finite array of point charges
%     arranged in the rocksalt structure (alternating F^- / Li^+ ions).
%     The potential energy V(x,y,z) = sum_i Q_probe * q_i / |r - r_i|
%     (Coulomb's law in atomic units) is evaluated on a dense grid above
%     the surface. The total potential is integrated over the active area
%     and converted from Hartree to eV.
%
%   Output:
%     potential_vs_height(z)  --  integrated potential at each height z
%
%   Dependencies:
%     calc_madelung_charge(row_extent, layer_count, probe_charge) returns the
%     3-D array of crystal point charges (see inline comment below).
%
%   Usage:
%     >> si_coulomb_potential
%     >> plot(height_vector, potential_vs_height);
%     >> xlabel('Height above surface (a.u.)');
%     >> ylabel('Integrated potential energy (eV)');
% =========================================================================

clear; clc;

%% ------------------------------------------------------------------------
%  CONFIGURATION
%  ------------------------------------------------------------------------

% --- Crystal lattice parameters ---

row_extent      = 2;        % Half-extent of crystal grid in x/y directions
                            % (total rows = 2*row_extent + 1)
layer_count     = 4;        % Number of ion layers below the surface
lattice_constant= 7.592;     % LiF rocksalt lattice constant (au, 2 × Ha = 2 × 3.7958)

% --- Probe parameters ---

point_charge    = -1;       % Charge of the probe particle (e.g. -1 for electron)

% --- Active surface scan grid ---

num_active_points = 500;    % Number of grid points per dimension on the
                            % active surface (total = N^2 points)
% Scan region: a square of side a_LiF/2 centred above the surface
active_site_x = -lattice_constant/4 : lattice_constant/num_active_points : lattice_constant/4;
active_site_y = -(-lattice_constant/4 : lattice_constant/num_active_points : lattice_constant/4);
% Area element for Riemann-sum integration over the active surface
area_element  = (lattice_constant / num_active_points)^2;

% --- Height (z) scan range ---

num_z_points  = 100;        % Number of z-sampling points
height_vector = linspace(2, 10, num_z_points);  % Heights 2..10 a.u.

%% ------------------------------------------------------------------------
%  BUILD CRYSTAL CHARGE DISTRIBUTION
%  ------------------------------------------------------------------------
%  Construct the discrete (2*row_extent+1) x (2*row_extent+1) x layer_count
%  grid of point charges. calc_madelung_charge returns a 3-D array containing the
%  alternating ionic charges (+1 for Li^+, -1 for F^-) at each site.
%  lattice_constant/2 is the nearest-neighbour ion spacing in rocksalt.

crystal_grid_x = (-row_extent:row_extent) * lattice_constant / 2;
crystal_grid_y = -1 * (-row_extent:row_extent) * lattice_constant / 2;
[coord_x, coord_y] = meshgrid(crystal_grid_x, crystal_grid_y);
charge_distribution = calc_madelung_charge(row_extent, layer_count, point_charge);

%% ------------------------------------------------------------------------
%  COMPUTE POTENTIAL ENERGY MAP AT EACH HEIGHT
%  ------------------------------------------------------------------------
%  For each probe height, we evaluate the Coulomb potential energy
%    V(ix,iy) = sum_{layers} sum_{crystal sites} Q * q_site / r
%  at every (active_site_x, active_site_y) position, then integrate over
%  the scan area to obtain the total potential energy for that height.

potential_vs_height = zeros(1, num_z_points);

for iZc = 1:num_z_points

    probe_height = height_vector(iZc);

    % 2-D potential energy map at current height z
    potential_energy_map = zeros(length(active_site_y), length(active_site_x));

    for ix = 1:length(active_site_x)
        for iy = 1:length(active_site_y)

            % Sum over crystal layers below the surface
            for nz = 1:layer_count
                % z-coordinate of the current layer (layer 1 = surface z=0,
                % subsequent layers at negative z below surface)
                coord_z = ones(2*row_extent+1, 2*row_extent+1) ...
                        * (1 - nz) * lattice_constant / 2;

                % Distance from probe to every crystal site in this layer
                r_distance = sqrt((coord_x - active_site_x(ix)).^2 ...
                                + (coord_y - active_site_y(iy)).^2 ...
                                + (coord_z - probe_height).^2);

                % Accumulate potential energy: V = sum_i Q_probe * q_i / |r - r_i|
                potential_energy_map(iy, ix) = potential_energy_map(iy, ix) ...
                    + sum(point_charge * charge_distribution(:, :, nz) ...
                          ./ r_distance, 'all');
            end
        end
    end

    % Integrate over active surface area and convert Hartree -> eV:
    %   total = [ sum_{x,y} V(x,y) * dS ] / (a/2)^2 * 27.2116
    % The division by (lattice_constant/2)^2 normalises per unit-cell area;
    % the factor 27.2116 eV/Ha converts from atomic energy units.
    total_potential = sum(potential_energy_map, 'all') ...
                    * area_element / (lattice_constant/2)^2 ...
                    * 27.2116;   % eV / Hartree

    potential_vs_height(iZc) = total_potential;

end

%% ------------------------------------------------------------------------
%  PLOT RESULT (optional — uncomment to visualise)
%  ------------------------------------------------------------------------

% plot(height_vector, potential_vs_height);
% xlabel('Height above surface (a.u.)');
% ylabel('Integrated potential energy (eV)');
% title('Coulomb potential energy vs. height above LiF surface');
