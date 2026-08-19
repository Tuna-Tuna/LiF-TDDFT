%% =============================================================================
%  calc_madelung_charge - Generate Madelung charge distribution for ionic crystal
%                         surface lattice
%
%  This function creates a 2D+layer Madelung charge array for an ionic crystal
%  lattice (e.g., LiF rock-salt structure). Alternating positive/negative charges
%  are arranged in a checkerboard pattern on each layer, with layer-parity sign
%  flip between adjacent layers.
%
%  USAGE:
%    charge_distribution = calc_madelung_charge(row_extent, layer_count, base_charge);
%
%  INPUTS:
%    row_extent   - Half-width of the 2D lattice plane (grid size = 2*row_extent+1)
%    layer_count  - Number of ion layers (z-direction)
%    base_charge  - Base charge magnitude (e.g., -1 for F⁻, +1 for Li⁺)
%
%  OUTPUTS:
%    charge_distribution - (2*row_extent+1) x (2*row_extent+1) x layer_count array
%                          with alternating ±base_charge values
%
%  REFERENCES:
%    Project: "Occupied-Space Constraints and Finite-Time Recovery in
%             Electron Detachment during F-/LiF(100) Scattering"
%    Code: Octopus TDDFT, real-time propagation
%
%  DEPENDENCIES: None (standalone function)
% =============================================================================

function charge_distribution = calc_madelung_charge(row_extent, layer_count, base_charge)

% Initialize template layer with all +base_charge
template_charge = ones(2*row_extent+1, 2*row_extent+1);

% Create alternating +/- pattern on first row of template
for i_column = 1:(2*row_extent+1)
    if mod(abs(i_column - row_extent - 1), 2) == 0
        template_charge(1, i_column) = 1;
    else
        template_charge(1, i_column) = -1;
    end
end

% Propagate the alternating pattern to all rows, flipping every other row
first_row_pattern = template_charge(1, :);
for i_row = 1:(2*row_extent+1)
    if mod(abs(row_extent - i_row - 1), 2) == 0
        template_charge(i_row, :) = first_row_pattern;
    else
        template_charge(i_row, :) = -1 * first_row_pattern;
    end
end

% Build 3D charge array with layer-parity sign flip
charge_distribution = ones(2*row_extent+1, 2*row_extent+1, layer_count);
for n_layer = 1:layer_count
    charge_distribution(:, :, n_layer) = template_charge * base_charge * (mod(n_layer, 2) * 2 - 1);
end

end
