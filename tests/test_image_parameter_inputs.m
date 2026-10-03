function test_image_parameter_inputs
% Interface validation only; the original image integration is not executed.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'matlab', 'analysis'));
try
    calc_image_potential(3, 0.2);
    error('test:unexpected', 'Missing dielectric inputs must fail.');
catch exception
    assert(strcmp(exception.identifier, 'imagepotential:parameters'));
end
try
    calc_image_potential(3, 0.2, NaN, 1);
    error('test:unexpected', 'Invalid dielectric inputs must fail.');
catch exception
    assert(~strcmp(exception.identifier, 'test:unexpected'));
end
fprintf('Image-potential parameter interface checks passed (2).\n');
end
