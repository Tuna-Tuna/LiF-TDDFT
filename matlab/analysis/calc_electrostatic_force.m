function result = calc_electrostatic_force(force_full_electronic, force_point_charge_reference)
%CALC_ELECTROSTATIC_FORCE Compatibility regression helper.
% The difference is a non-unique electronic contribution, not a unique
% Pauli/exchange force decomposition.
if ~isequal(size(force_full_electronic), size(force_point_charge_reference))
    error('Force arrays must have the same shape.');
end
result.force_full_electronic = force_full_electronic;
result.force_point_charge_reference = force_point_charge_reference;
result.force_full_minus_pc = force_full_electronic - force_point_charge_reference;
end
