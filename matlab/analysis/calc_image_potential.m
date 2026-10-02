function image_potential = calc_image_potential(z_grid, velocity_list)
% CALC_IMAGE_POTENTIAL  Compute image potential energy for a charged
%                       projectile moving near a dielectric surface.
%
%   Computes the dynamical image potential experienced by a moving point
%   charge above a LiF(001) surface using the specular-reflection model
%   (Echenique-Pendry formulation). The potential arises from the
%   polarization response of the dielectric, which induces a time-dependent
%   image charge that lags behind the projectile.
%
%   PHYSICAL MODEL:
%     The image potential V_image(z, v) is computed via numerical
%     integration of the dielectric response function over momentum (W) and
%     frequency (T) space. The key quantities are:
%
%       E1(W) = eps_d * eps_0 / (eps_0^2 - W^2)
%       E2(W, T) = (1/(T^2-1)) * eps_d * eps_0 / (eps_0^2 - (T*W)^2)
%       Re(W) = (E1^2 + E2^2 - 1) / ((E1+1)^2 + E2^2)   (reflection coeff)
%       V(z,v) = -1/(pi*v) * integral( K0(2*W*z/v) * Re(W) ) dW
%
%     where eps_0 = 17.1 eV and eps_d = 14.9 eV are characteristic
%     energies of the LiF dielectric function.
%
%   INPUT:
%       z_grid       - Height grid [num_z x 1] in atomic units (Bohr).
%       velocity_list - Projectile velocity list [1 x num_v] in atomic units.
%
%   OUTPUT:
%       image_potential - Matrix [num_z x num_v] of image potential energy
%                         values in eV. Column j corresponds to velocity
%                         velocity_list(j), row i to height z_grid(i).
%
%   USAGE:
%       z = linspace(1, 10, 1801);
%       v = [0.1, 0.15, 0.2, 0.3, 0.4, 0.5];
%       V_image = calc_image_potential(z, v);
%       [~, F_image] = gradient(V_image, 0.005);
%       F_image = -F_image / 27.2116;  % Convert to atomic force units
%
%   REFERENCES:
%     - Echenique, P.M. & Pendry, J.B. J. Phys. C 8, 2936 (1975)
%
%   DEPENDENCIES:
%     - MATLAB (besselk for modified Bessel function K0)
%
%   NOTE:
%     The dielectric parameters (E0 = 17.1 eV, Ed = 14.9 eV) are fitted
%     to the LiF optical data. The integration cutoffs (W_max = 0.5,
%     T_max = 0.5) are chosen based on convergence testing for the
%     relevant velocity range (v ~ 0.1-0.5 au).
%
%   See also: besselk, trapz, gradient

%% ==================== Physical Constants ====================
HARTREE_TO_EV    = 27.2116;          % eV / Hartree
DIELEC_ENERGY_D  = 14.9 / HARTREE_TO_EV;  % Ed, Hartree
DIELEC_ENERGY_0  = 17.1 / HARTREE_TO_EV;  % E0, Hartree

%% ==================== Integration Parameters ====================
% Momentum integration grid (W-space)
W_min   = 1e-5;
W_max   = 0.5;
W_step  = 1e-5;
W_grid  = (W_min : W_step : W_max)';   % Momentum values in Hartree

% Frequency integration grid (T = omega/omega_s, dimensionless)
T_min   = 1e-5;
T_max   = 0.5;
T_step  = 1e-5;
T_grid  = (T_min : T_step : T_max)';   % Dimensionless frequency ratio

%% ==================== Pre-compute Dielectric Response ====================
% E1: Real part of the surface dielectric function (energy-dependent)
%     E1(W) = Ed * E0 / (E0^2 - W^2)
dielectric_real = DIELEC_ENERGY_D * DIELEC_ENERGY_0 ./ ...
    (DIELEC_ENERGY_0^2 - W_grid.^2);

% E2 coefficient (W-dependent prefactor for T-integral)
%     Numerator of the imaginary part integrand
E2_numerator = DIELEC_ENERGY_D * DIELEC_ENERGY_0;

num_z = length(z_grid);
num_v = length(velocity_list);

image_potential_hartree = zeros(num_z, num_v);
image_potential = zeros(num_z, num_v);

%% ==================== Compute Image Potential ====================
for idx_v = 1:num_v
    projectile_velocity = velocity_list(idx_v);

    for idx_z = 1:num_z
        probe_height = z_grid(idx_z);

        % ---- Frequency integration (T-space): E2 component ----
        % E2(W, T) = 1/(T^2-1) * (Ed * E0) / (E0^2 - (T*W)^2)
        % This is a 2D integral over both W and T.

        % Build the integrand over T for each W point, then integrate
        % over T via trapz
        T_difference_sq = T_grid.^2 - 1;                  % T^2 - 1
        denominator_T = DIELEC_ENERGY_0^2 - (T_grid .* W_grid').^2;  % E0^2 - (T*W)^2
        integrand_T = (1 ./ T_difference_sq') .* (E2_numerator ./ denominator_T);

        dielectric_imag = -2 / (pi * projectile_velocity) .* ...
            trapz(T_grid, integrand_T, 1);               % Integrate over T (dim 1)

        % ---- Reflection coefficient ----
        % Re(W) = (E1^2 + E2^2 - 1) / ((E1 + 1)^2 + E2^2)
        numerator_reflection = dielectric_real.^2 + dielectric_imag.^2 - 1;
        denominator_reflection = (dielectric_real + 1).^2 + dielectric_imag.^2;
        reflection_coeff = numerator_reflection ./ denominator_reflection;

        % ---- Momentum integration (W-space): Kernel and final integral ----
        % K0(2 * W * z / v) is the modified Bessel function of second kind
        bessel_argument = 2 * W_grid .* (probe_height / projectile_velocity);
        bessel_kernel = besselk(0, bessel_argument);

        % The full integrand
        integrand_W = bessel_kernel .* reflection_coeff;

        % Trapz integration over W
        integral_result = trapz(W_grid, integrand_W);

        % Final potential in Hartree
        image_potential_hartree(idx_z, idx_v) = -1 / (pi * projectile_velocity) ...
            * integral_result;

        % Convert to eV
        image_potential(idx_z, idx_v) = image_potential_hartree(idx_z, idx_v) ...
            * HARTREE_TO_EV;
    end
end

%% ==================== Post-processing Note ====================
% To obtain the image force in atomic units, compute:
%   [~, F_image] = gradient(image_potential, z_step);
%   F_image = -F_image / HARTREE_TO_EV;
%
% This is typically done in the calling script (run_trajectory_simulation.m).

end
