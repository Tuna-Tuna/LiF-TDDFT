function [density_first, density_second, cell_volume] = compare_on_native_grids(first, second, first_axes, second_axes)
% Compare stored (z,y,x) densities on matching uniform (x,y,z) coordinate axes.
% No intermediate samples or zero-filled missing densities are constructed.
spacing = zeros(1, 3);
if numel(first_axes) ~= 3 || numel(second_axes) ~= 3
    error('nativegrid:axes', 'Three recorded coordinate axes are required.');
end
for dim = 1:3
    a = first_axes{dim}(:);
    b = second_axes{dim}(:);
    if numel(a) < 2 || numel(a) ~= numel(b) || ...
            any(~isfinite(a)) || any(~isfinite(b)) || any(abs(a-b) > 1e-10)
        error('nativegrid:matching', 'Grids must contain matching recorded coordinates.');
    end
    increments = diff(a);
    if any(increments <= 0) || any(abs(increments-increments(1)) > 1e-10)
        error('nativegrid:uniform', 'Density integration requires a recorded uniform grid.');
    end
    spacing(dim) = increments(1);
end
expected_size = [numel(first_axes{3}), numel(first_axes{2}), numel(first_axes{1})];
if ~isequal(size(first), expected_size) || ~isequal(size(second), expected_size)
    error('nativegrid:shape', 'Density dimensions must agree with the recorded axes.');
end
if any(~isfinite(first(:))) || any(~isfinite(second(:)))
    error('nativegrid:missing', 'Source density contains missing values; zero filling is forbidden.');
end
cell_volume = prod(spacing);
density_first = permute(first, [2, 3, 1]);
density_second = permute(second, [2, 3, 1]);
end
