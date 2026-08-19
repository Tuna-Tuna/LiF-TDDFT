% octave_startup.m - GNU Octave compatibility setup
%
% Run this once at the start of your Octave session before using any
% analysis scripts that read NetCDF density files:
%
%   run('matlab/utils/octave_startup.m');
%
% This script loads required Octave packages and verifies compatibility.
%
% NOTE: Octave >= 4.4 is required for full compatibility (contains(),
%       griddedInterpolant with cell arrays, etc.). Octave 6.0+ recommended.

fprintf('=== GNU Octave Compatibility Setup ===\n');

% Load netcdf package (required for ncread in analysis scripts)
try
    pkg load netcdf;
    fprintf('  netcdf package loaded successfully.\n');
catch
    warning('netcdf package not found. Install it with:');
    warning('  pkg install -forge netcdf');
    warning('  Then re-run: pkg load netcdf');
end

fprintf('=== Setup Complete ===\n');
