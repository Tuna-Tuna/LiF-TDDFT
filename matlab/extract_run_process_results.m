function manifest = extract_run_process_results(inputMat, outputDir)
%EXTRACT_RUN_PROCESS_RESULTS Export the saved run_process workspace safely.
%   This function reads only named variables from PdE5.mat. It deliberately
%   excludes allCoordinates (~760 MB), interpolant objects, local paths, and
%   the last-loaded data struct. Numeric arrays are exported as CSV so the
%   Python revision workflow does not require scipy.io.

if nargin < 1 || strlength(inputMat) == 0
    inputMat = "PdE5.mat";
end
if nargin < 2 || strlength(outputDir) == 0
    outputDir = fullfile("data", "processed", "run_process");
end
inputMat = char(inputMat);
outputDir = char(outputDir);
if ~isfile(inputMat)
    error('Input MAT file not found: %s', inputMat);
end
if ~isfolder(outputDir)
    mkdir(outputDir);
end

names = {
    'vcal', 'zcal', 'result_hsemt', 'result_esemt', 'result_loss', ...
    'ZImP', 'Image', 'Fimage', 'z', 'zdeltaE_Cap', ...
    'Zhr', 'Zhrx', 'zPcapture', 'zPloss', 'zdE', 'Pfin', ...
    'last_non_zero', 'Zturn', 'E', 'expe', 'coor_ac_x', 'force_vh', ...
    'a_LiF', 'alphaPositive', 'alphaNegtive', 'MadField_fix', ...
    'dEF', 'gama', 'FM', 'alpha', 'angstrom'
};
available = whos('-file', inputMat);
availableNames = {available.name};
missing = setdiff(names, availableNames);
if ~isempty(missing)
    error('Required variables missing from MAT file: %s', strjoin(missing, ', '));
end

loaded = load(inputMat, names{:});
arrayNames = {
    'vcal', 'zcal', 'result_hsemt', 'result_esemt', 'result_loss', ...
    'ZImP', 'Image', 'Fimage', 'z', 'zdeltaE_Cap', ...
    'Zhr', 'Zhrx', 'zPcapture', 'zPloss', 'zdE', 'Pfin', ...
    'last_non_zero', 'Zturn', 'E', 'expe', 'coor_ac_x', 'force_vh'
};
records = struct('name', {}, 'shape', {}, 'class', {}, 'file', {});
for k = 1:numel(arrayNames)
    name = arrayNames{k};
    value = loaded.(name);
    fileName = [name '.csv'];
    writematrix(value, fullfile(outputDir, fileName));
    records(k).name = name;
    records(k).shape = size(value);
    records(k).class = class(value);
    records(k).file = fileName;
end

scalarNames = {
    'a_LiF', 'alphaPositive', 'alphaNegtive', 'MadField_fix', ...
    'dEF', 'gama', 'FM', 'alpha', 'angstrom'
};
parameters = struct();
for k = 1:numel(scalarNames)
    parameters.(scalarNames{k}) = loaded.(scalarNames{k});
end

manifest = struct();
manifest.schema_version = 1;
manifest.source_mat = char(java.io.File(inputMat).getCanonicalPath());
sourceInfo = dir(inputMat);
manifest.source_bytes = sourceInfo.bytes;
manifest.excluded_variables = {
    'allCoordinates', 'data', 'imageInterpolants', 'fimageInterpolants', ...
    'forceInterpolants', 'result_lossInterpolants', 'cap_data', 'los_data'
};
manifest.arrays = records;
manifest.parameters = parameters;
manifest.notes = {
    'zPloss and zPcapture are saved legacy per-event probabilities.', ...
    'result_loss may include values outside [0,1]; it is exported unchanged and Python does not clip or rescale it.', ...
    'No missing values or out-of-domain events may be silently extrapolated.'
};
jsonText = jsonencode(manifest, PrettyPrint=true);
fid = fopen(fullfile(outputDir, 'manifest.json'), 'w', 'n', 'UTF-8');
if fid < 0
    error('Cannot write extraction manifest.');
end
cleanup = onCleanup(@() fclose(fid));
fwrite(fid, jsonText, 'char');
fwrite(fid, newline, 'char');
fprintf('Exported %d arrays to %s\n', numel(arrayNames), outputDir);
end
