function manifest = extract_classical_force_fields(inputMat, outputDir)
%EXTRACT_CLASSICAL_FORCE_FIELDS Export compact force fields from PdE5.mat.
%   The saved allCoordinates cell array is large, but the legacy classical
%   trajectory loop only consumes one projectile-force column and one path
%   coordinate per height and velocity.  This function exports that compact
%   subset without copying the full coordinate matrices into the Python
%   workflow.

if nargin < 1 || strlength(inputMat) == 0
    inputMat = "PdE5.mat";
end
if nargin < 2 || strlength(outputDir) == 0
    outputDir = fullfile("data", "processed", "run_process", "force_fields");
end
inputMat = char(inputMat);
outputDir = char(outputDir);
if ~isfile(inputMat)
    error('Input MAT file not found: %s', inputMat);
end
if ~isfolder(outputDir)
    mkdir(outputDir);
end

available = whos('-file', inputMat);
availableNames = {available.name};
required = {'allCoordinates', 'vcal', 'zcal'};
missing = setdiff(required, availableNames);
if ~isempty(missing)
    error('Required variables missing from MAT file: %s', strjoin(missing, ', '));
end

saved = load(inputMat, required{:});
allCoordinates = saved.allCoordinates;
vcal = saved.vcal(:)';
zcal = saved.zcal(:)';
nHeights = size(allCoordinates, 1);
nVelocities = size(allCoordinates, 2);
nSurfaceIons = 18;
if nHeights ~= numel(zcal) || nVelocities ~= numel(vcal)
    error('allCoordinates dimensions do not match zcal/vcal.');
end

records = struct('velocity_au', {}, 'rows', {}, 'ps', {}, 'pe', {}, 'file', {});
for velocityIndex = 1:nVelocities
    coordinates = allCoordinates{1, velocityIndex};
    if isempty(coordinates)
        error('Missing allCoordinates{%d,%d}.', 1, velocityIndex);
    end
    ps = find(coordinates(:, 3) - coordinates(:, 6) < 0, 1, 'first');
    pe = find(coordinates(:, 3) - coordinates(:, 6 + 3*6) < 0, 1, 'first');
    if isempty(ps) || isempty(pe) || pe < ps
        error('Could not determine legacy path window for velocity index %d.', velocityIndex);
    end

    forceByHeight = zeros(pe - ps + 1, nHeights);
    pathX = [];
    for heightIndex = 1:nHeights
        coordinates = allCoordinates{heightIndex, velocityIndex};
        if isempty(coordinates) || size(coordinates, 1) < pe
            error('Incomplete coordinates at height %d, velocity %d.', heightIndex, velocityIndex);
        end
        forceByHeight(:, heightIndex) = coordinates(ps:pe, end - 3*nSurfaceIons);
        if heightIndex == nHeights
            pathX = coordinates(ps:pe, 3);
        end
    end

    variableNames = [{'path_x_bohr'}, compose('force_z_h%g_eV_A', zcal)];
    tableData = array2table([pathX, forceByHeight], 'VariableNames', variableNames);
    velocityLabel = strrep(sprintf('%.6g', vcal(velocityIndex)), '.', 'p');
    fileName = sprintf('force_field_v%s.csv', velocityLabel);
    writetable(tableData, fullfile(outputDir, fileName));

    records(velocityIndex).velocity_au = vcal(velocityIndex);
    records(velocityIndex).rows = size(forceByHeight, 1);
    records(velocityIndex).ps = ps;
    records(velocityIndex).pe = pe;
    records(velocityIndex).file = fileName;
end

manifest = struct();
manifest.schema_version = 1;
manifest.source_mat = char(java.io.File(inputMat).getCanonicalPath());
sourceInfo = dir(inputMat);
manifest.source_bytes = sourceInfo.bytes;
manifest.zcal_bohr = zcal;
manifest.projectile_force_column_rule = 'end - 3*Npar, Npar=18';
manifest.path_window_rule = 'legacy processDataForVelocity in run_process.m';
manifest.velocity_force_fields = records;
manifest.notes = {
    'Files contain only the force column and path coordinate used by the legacy classical loop.', ...
    'The full allCoordinates matrices are not exported.', ...
    'These derived force fields are intended for deterministic replay and audit, not new Octopus production.'
};

jsonText = jsonencode(manifest, PrettyPrint=true);
fid = fopen(fullfile(outputDir, 'manifest.json'), 'w', 'n', 'UTF-8');
if fid < 0
    error('Cannot write force-field extraction manifest.');
end
cleanup = onCleanup(@() fclose(fid));
fwrite(fid, jsonText, 'char');
fwrite(fid, newline, 'char');
fprintf('Exported %d compact force fields to %s\n', nVelocities, outputDir);
end
