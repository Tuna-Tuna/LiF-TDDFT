function summary = run_nonproduction_p0_analysis(repoRoot, archiveRoot)
%RUN_NONPRODUCTION_P0_ANALYSIS Execute the revision P0 analysis only.
%   This routine consumes archived result tables and PdE5-derived arrays.
%   It does not invoke Octopus or alter any production calculation.

if nargin < 1 || strlength(repoRoot) == 0
    repoRoot = "F:/github/LiF-TDDFT";
end
if nargin < 2 || strlength(archiveRoot) == 0
    archiveRoot = "F:/codex/JCTC_纯分析结果与支撑数据_20260727";
end
repoRoot = char(repoRoot);
archiveRoot = char(archiveRoot);

processedRoot = fullfile(repoRoot, 'data', 'processed', 'run_process');
forceFieldRoot = fullfile(processedRoot, 'force_fields');
tableRoot = fullfile(repoRoot, 'results', 'tables');
figureRoot = fullfile(repoRoot, 'results', 'figures');
if ~isfolder(tableRoot); mkdir(tableRoot); end
if ~isfolder(figureRoot); mkdir(figureRoot); end

vcal = readmatrix(fullfile(processedRoot, 'vcal.csv'));
vcal = vcal(:)';
zcal = readmatrix(fullfile(processedRoot, 'zcal.csv'));
zcal = zcal(:)';
ZhrSaved = readmatrix(fullfile(processedRoot, 'Zhr.csv'));
ZturnSaved = readmatrix(fullfile(processedRoot, 'Zturn.csv'));
zPloss = readmatrix(fullfile(processedRoot, 'zPloss.csv'));
zPcaptureSaved = readmatrix(fullfile(processedRoot, 'zPcapture.csv'));
zdE = readmatrix(fullfile(processedRoot, 'zdE.csv'));
PfinSaved = readmatrix(fullfile(processedRoot, 'Pfin.csv'));
manifest = jsondecode(fileread(fullfile(processedRoot, 'manifest.json')));
params = manifest.parameters;

%% 1. Reconstruct the legacy classical trajectories with normal velocity.
eventRows = cell(nnz(ZhrSaved ~= 0 & isfinite(ZhrSaved)), 14);
eventCursor = 0;
replayRows = cell(numel(vcal), 13);
for iv = 1:numel(vcal)
    label = strrep(sprintf('%.6g', vcal(iv)), '.', 'p');
    forcePath = fullfile(forceFieldRoot, sprintf('force_field_v%s.csv', label));
    forceTable = readtable(forcePath, 'VariableNamingRule', 'preserve');
    pathX = forceTable{:, 1};
    forceByHeight = forceTable{:, 2:end};
    splines = cell(size(forceByHeight, 1), 1);
    for row = 1:size(forceByHeight, 1)
        splines{row} = griddedInterpolant(zcal, forceByHeight(row, :), 'spline');
    end

    alpha = params.alpha;
    angstrom = params.angstrom;
    mass = params.FM;
    vx = -vcal(iv) * cos(alpha);
    vz = -vcal(iv) * sin(alpha);
    deltaT = diff(pathX) * angstrom / vx;
    currentHeight = 10.0;
    minimumHeight = currentHeight;
    legacyMinimumHeight = currentHeight;
    eventIndex = 0;
    safetyCycles = 0;
    reconstructedHeights = [];
    reconstructedNormalVelocity = [];
    midpoint = round(numel(pathX) / 2);
    while isreal(currentHeight) && currentHeight <= 10.0
        safetyCycles = safetyCycles + 1;
        if safetyCycles > 1000
            error('Trajectory replay exceeded safety limit for v=%g.', vcal(iv));
        end
        legacyMinimumHeight = min(legacyMinimumHeight, currentHeight);
        for izi = 2:numel(pathX)
            forceAu = splines{izi}(currentHeight) / (27.2116 * angstrom);
            acceleration = forceAu / mass;
            dt = deltaT(izi - 1);
            currentHeight = currentHeight + vz * dt + 0.5 * acceleration * dt^2;
            vz = vz + acceleration * dt;
            minimumHeight = min(minimumHeight, currentHeight);
            if izi == midpoint
                eventIndex = eventIndex + 1;
                reconstructedHeights(eventIndex, 1) = currentHeight; %#ok<AGROW>
                reconstructedNormalVelocity(eventIndex, 1) = vz; %#ok<AGROW>
                withinDatabase = currentHeight >= min(zcal) && currentHeight <= max(zcal);
                belowDatabase = currentHeight < min(zcal);
                aboveDatabase = currentHeight > max(zcal);
                eventCursor = eventCursor + 1;
                eventRows(eventCursor, :) = { ...
                    sprintf('legacy_v%.2f', vcal(iv)), 1.0, iv, eventIndex, ...
                    vcal(iv), currentHeight, abs(vx), vz, abs(vz)/abs(vx), ...
                    0.0, withinDatabase, belowDatabase, aboveDatabase, ...
                    (eventIndex - 0.5) * params.a_LiF};
            end
        end
    end

    savedHeights = ZhrSaved(iv, :);
    savedHeights = savedHeights(savedHeights ~= 0 & isfinite(savedHeights));
    countMatch = numel(savedHeights) == numel(reconstructedHeights);
    if countMatch
        maxHeightError = max(abs(savedHeights(:) - reconstructedHeights(:)));
    else
        maxHeightError = NaN;
    end
    savedFinalIndex = find(PfinSaved(iv, :) ~= 0, 1, 'last');
    replayRows(iv, :) = {vcal(iv), numel(reconstructedHeights), numel(savedHeights), ...
        countMatch, maxHeightError, minimumHeight, ZturnSaved(iv), ...
        legacyMinimumHeight, abs(legacyMinimumHeight - ZturnSaved(iv)), ...
        minimumHeight - ZturnSaved(iv), max(abs(reconstructedNormalVelocity)), ...
        max(abs(reconstructedNormalVelocity))/abs(vx), savedFinalIndex};
end
if eventCursor ~= size(eventRows, 1)
    error('Reconstructed event count does not match the saved nonzero event count.');
end

events = cell2table(eventRows, 'VariableNames', { ...
    'trajectory_id', 'trajectory_weight', 'velocity_index', 'event_index', ...
    'incident_velocity_au', 'surface_height_bohr', 'v_parallel_au', ...
    'v_perpendicular_au', 'abs_vperp_over_vparallel', ...
    'minimum_lateral_distance_bohr', 'within_td_database', ...
    'below_td_database', 'above_td_database', 'path_coordinate_bohr'});
writetable(events, fullfile(tableRoot, 'classical_event_reconstruction.csv'));

replay = cell2table(replayRows, 'VariableNames', { ...
    'incident_velocity_au', 'reconstructed_event_count', 'saved_event_count', ...
    'event_count_match', 'max_event_height_error_bohr', ...
    'true_substep_turning_height_bohr', 'saved_turning_height_bohr', ...
    'reconstructed_legacy_turning_height_bohr', 'legacy_turning_height_absolute_error_bohr', ...
    'true_minus_saved_turning_height_bohr', 'max_abs_vperpendicular_au', ...
    'max_abs_vperpendicular_over_vparallel', 'saved_final_probability_index'});
writetable(replay, fullfile(tableRoot, 'classical_replay_validation.csv'));

thresholds = [3.5, 3.0, 2.5, 2.0];
wgeomRows = cell(numel(vcal) + 1, 12);
for row = 1:(numel(vcal) + 1)
    if row <= numel(vcal)
        mask = events.velocity_index == row;
        velocityValue = vcal(row);
        scope = sprintf('v=%.2f', velocityValue);
        incidentWeight = 1.0;
    else
        mask = true(height(events), 1);
        velocityValue = NaN;
        scope = 'equal-weight_six-velocity_ensemble';
        incidentWeight = numel(vcal);
    end
    selected = events(mask, :);
    weights = selected.trajectory_weight;
    meanCollisions = sum(weights) / incidentWeight;
    belowFractions = arrayfun(@(h) sum(weights(selected.surface_height_bohr < h))/incidentWeight, thresholds);
    covered = sum(weights(selected.within_td_database));
    total = sum(weights);
    wgeomRows(row, :) = {scope, velocityValue, height(selected), meanCollisions, ...
        belowFractions(1), belowFractions(2), belowFractions(3), belowFractions(4), ...
        covered, total-covered, covered/total, ...
        max(selected.abs_vperp_over_vparallel)};
end
wgeomSummary = cell2table(wgeomRows, 'VariableNames', { ...
    'scope', 'incident_velocity_au', 'event_count', 'mean_effective_collision_count', ...
    'events_per_incident_below_h3p5', 'events_per_incident_below_h3p0', ...
    'events_per_incident_below_h2p5', 'events_per_incident_below_h2p0', ...
    'covered_event_weight', 'out_of_domain_event_weight', ...
    'td_database_coverage_fraction', 'max_abs_vperp_over_vparallel'});
writetable(wgeomSummary, fullfile(tableRoot, 'wgeom_diagnostic_summary.csv'));

heightEdges = [1, 1.5, 2, 2.5, 3, 3.5, 4, 5, 7, 10.1];
histRows = cell(numel(vcal) * (numel(heightEdges) - 1), 6);
histCursor = 0;
for iv = 1:numel(vcal)
    mask = events.velocity_index == iv;
    counts = histcounts(events.surface_height_bohr(mask), heightEdges);
    for ib = 1:numel(counts)
        histCursor = histCursor + 1;
        histRows(histCursor, :) = {vcal(iv), heightEdges(ib), heightEdges(ib+1), ...
            counts(ib), counts(ib)/sum(counts), counts(ib)};
    end
end
wgeomBins = cell2table(histRows, 'VariableNames', { ...
    'incident_velocity_au', 'height_bin_lower_bohr', 'height_bin_upper_bohr', ...
    'event_count', 'within_velocity_event_fraction', 'events_per_incident'});
writetable(wgeomBins, fullfile(tableRoot, 'wgeom_height_bins.csv'));

%% 2. Align continuous Q20 descriptors with the v=0.30 TD trajectory.
q20File = findOne(archiveRoot, 'Q20_continuous_metrics_20260731.csv');
q20 = readtable(q20File, 'VariableNamingRule', 'preserve');
abFile = findOne(archiveRoot, 'h35_v030_spatial_raw.csv');
ab = readtable(abFile, 'VariableNamingRule', 'preserve');
backgroundFile = findOne(fullfile(archiveRoot, '04_背景参考与脱附计数'), 'v030_spatial_raw.csv');
background = readtable(backgroundFile, 'VariableNamingRule', 'preserve');
forceFile = findOne(archiveRoot, 'onlyForces', fullfile('scalar_sources', 'v030'));
coordinateFile = findOne(archiveRoot, 'onlyCoordinates', fullfile('scalar_sources', 'v030'));
forces = readOctopusNumeric(forceFile);
coordinates = readOctopusNumeric(coordinateFile);
if size(forces, 1) ~= size(coordinates, 1) || any(forces(:,1) ~= coordinates(:,1))
    error('v030 force and coordinate series are not aligned.');
end
projectileX = coordinates(:, 3);
projectileForce = forces(:, 3:5);

qX = q20.projectile_x_A;
alignedForceX = interpDescending(projectileX, projectileForce(:,1), qX);
alignedForceY = interpDescending(projectileX, projectileForce(:,2), qX);
alignedForceZ = interpDescending(projectileX, projectileForce(:,3), qX);
alignedForceMagnitude = sqrt(alignedForceX.^2 + alignedForceY.^2 + alignedForceZ.^2);

abNormLoss = ab.norm(1) - interpDescending(ab.x, ab.norm, qX);
bgNormLoss = background.norm(1) - interpDescending(background.x, background.norm, qX);
interactionNormLoss = abNormLoss - bgNormLoss;
abCapGrowth = interpDescending(ab.x, ab.cap_population, qX) - ab.cap_population(1);
bgCapGrowth = interpDescending(background.x, background.cap_population, qX) - background.cap_population(1);
interactionCapGrowth = abCapGrowth - bgCapGrowth;
abR35Loss = ab.('pop_r3.5')(1) - interpDescending(ab.x, ab.('pop_r3.5'), qX);
bgR35Loss = background.('pop_r3.5')(1) - interpDescending(background.x, background.('pop_r3.5'), qX);
interactionR35Loss = abR35Loss - bgR35Loss;

q20Td = table(q20.point_id, q20.RFF_A, q20.projectile_x_A, q20.projectile_z_A, ...
    q20.overlap_max, q20.overlap_frobenius, q20.lowdin_half_L1_e, ...
    q20.vW_delta_cutoff_1e10_eV, alignedForceX, alignedForceY, alignedForceZ, ...
    alignedForceMagnitude, abNormLoss, bgNormLoss, interactionNormLoss, ...
    interactionCapGrowth, interactionR35Loss, ...
    'VariableNames', {'point_id','RFF_A','projectile_x_A','projectile_z_A', ...
    'overlap_max','overlap_frobenius','lowdin_half_L1_e','vW_delta_eV', ...
    'td_force_x_eV_A','td_force_y_eV_A','td_force_z_eV_A','td_force_magnitude_eV_A', ...
    'ab_norm_loss_e','background_norm_loss_e','interaction_induced_norm_loss_e', ...
    'interaction_induced_cap_growth_e','interaction_induced_r3p5_population_loss_e'});
writetable(q20Td, fullfile(tableRoot, 'q20_td_common_RFF_alignment.csv'));

q20Correlation = struct();
q20Correlation.points = height(q20Td);
q20Correlation.log_overlap_vs_force_pearson = pearson(log10(q20Td.overlap_frobenius), q20Td.td_force_magnitude_eV_A);
q20Correlation.log_lowdin_rearrangement_vs_force_pearson = pearson(log10(q20Td.lowdin_half_L1_e), q20Td.td_force_magnitude_eV_A);
q20Correlation.log_vw_proxy_vs_force_pearson = pearson(log10(q20Td.vW_delta_eV), q20Td.td_force_magnitude_eV_A);
q20Correlation.force_peak_point = char(q20Td.point_id{end});
q20Correlation.force_peak_magnitude_eV_A = q20Td.td_force_magnitude_eV_A(end);
writeJson(fullfile(tableRoot, 'q20_td_alignment_summary.json'), q20Correlation);

figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1300 900]);
tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
nexttile;
semilogy(q20Td.RFF_A, q20Td.overlap_frobenius, 'o-', 'LineWidth', 1.5); hold on;
semilogy(q20Td.RFF_A, q20Td.lowdin_half_L1_e, 's-', 'LineWidth', 1.5);
set(gca, 'XDir', 'reverse'); grid on;
xlabel('R_{F-F} (Å)'); ylabel('Orthogonalization descriptor');
legend('||S_{AB}||_F', 'Löwdin half-L1 (e)', 'Location', 'best');
title('Continuous Q20 descriptors');
nexttile;
plot(q20Td.RFF_A, q20Td.vW_delta_eV, 'o-', 'LineWidth', 1.5); hold on;
plot(q20Td.RFF_A, q20Td.td_force_magnitude_eV_A, 's-', 'LineWidth', 1.5);
set(gca, 'XDir', 'reverse'); grid on;
xlabel('R_{F-F} (Å)'); ylabel('eV or eV/Å');
legend('\DeltaT_{vW} proxy', '|F_{TD}|', 'Location', 'best');
title('Proxy and TD force on common R');
nexttile;
plot(q20Td.RFF_A, q20Td.interaction_induced_norm_loss_e, 'o-', 'LineWidth', 1.5); hold on;
plot(q20Td.RFF_A, q20Td.interaction_induced_cap_growth_e, 's-', 'LineWidth', 1.5);
set(gca, 'XDir', 'reverse'); grid on;
xlabel('R_{F-F} (Å)'); ylabel('Interaction-induced change (e)');
legend('KS norm loss', 'CAP population growth', 'Location', 'best');
title('Incoming-branch electronic loss diagnostics');
nexttile;
bar(categorical(string(wgeomSummary.scope(1:6))), wgeomSummary.td_database_coverage_fraction(1:6));
ylim([0 1.05]); grid on; ylabel('TD database coverage'); xtickangle(35);
title('Legacy event coverage by velocity');
exportgraphics(gcf, fullfile(figureRoot, 'p0_q20_td_wgeom_summary.png'), 'Resolution', 240);
close(gcf);

%% 3. Quantify the historical near-constant trajectory protocol.
protocolSpeeds = [0.20, 0.30, 0.40];
protocolRows = cell(numel(protocolSpeeds), 12);
for iv = 1:numel(protocolSpeeds)
    speedCode = sprintf('v%03d', round(protocolSpeeds(iv)*100));
    spatialFile = findOne(archiveRoot, sprintf('h35_%s_spatial_raw.csv', speedCode));
    spatial = readtable(spatialFile, 'VariableNamingRule', 'preserve');
    velocityConversion = 27.2116 / params.angstrom;
    targetVxStored = -protocolSpeeds(iv) * velocityConversion;
    productionStart = find(abs(spatial.vx - targetVxStored)/abs(targetVxStored) <= 0.01, 1, 'first');
    if isempty(productionStart)
        error('Could not identify production segment for %s.', speedCode);
    end
    production = productionStart:height(spatial);
    z0 = spatial.z(productionStart);
    protocolRows(iv, :) = {speedCode, protocolSpeeds(iv), productionStart, ...
        spatial.iteration(productionStart), targetVxStored/velocityConversion, ...
        spatial.vx(productionStart)/velocityConversion, spatial.vx(end)/velocityConversion, ...
        max(abs(spatial.vx(production) - targetVxStored))/abs(targetVxStored), ...
        max(abs(spatial.vz(production)))/velocityConversion, ...
        max(abs(spatial.z(production) - z0)), spatial.x(productionStart), spatial.x(end)};
end
protocol = cell2table(protocolRows, 'VariableNames', { ...
    'case_id','nominal_velocity_au','production_start_row','production_start_iteration', ...
    'target_vx_au','production_start_vx_au','final_vx_au', ...
    'max_relative_vx_deviation_from_target','max_abs_vz_au','max_abs_height_drift_A', ...
    'production_start_x_A','final_x_A'});
writetable(protocol, fullfile(tableRoot, 'historical_td_protocol_diagnostics.csv'));

%% 4. SI Eq. (S14) audit and parameter sensitivity on saved encounter geometry.
% The sensitivity loop re-evaluates S14 at alternate physical parameters.
% It never rescales, normalizes, or replaces the baseline probability data.
gammaMultipliers = [0.8 0.9 1.0 1.1 1.2];
energyShifts = [-0.02 -0.01 0 0.01 0.02];
sensitivityRows = cell(numel(vcal) * (numel(gammaMultipliers) + numel(energyShifts)), 7);
sensitivityCursor = 0;
baselineCaptureMaxError = 0;
for iv = 1:numel(vcal)
    valid = find(ZhrSaved(iv,:) ~= 0 & isfinite(ZhrSaved(iv,:)));
    heights = ZhrSaved(iv, valid); %#ok<NASGU> retained for event audit output
    detachment = zPloss(iv, valid);
    if any(~isfinite(detachment), 'all') || any(detachment < 0 | detachment > 1, 'all')
        error('Saved detachment probability is invalid; clipping or rescaling is forbidden.');
    end
    baseEnergy = zdE(iv, valid);
    baseCapture = demkovS14(vcal(iv), baseEnergy, params.gama);
    baselineCaptureMaxError = max(baselineCaptureMaxError, max(abs(baseCapture - zPcaptureSaved(iv, valid))));
    for multiplier = gammaMultipliers
        capture = demkovS14(vcal(iv), baseEnergy, params.gama * multiplier);
        finalYield = propagate(detachment, capture);
        sensitivityCursor = sensitivityCursor + 1;
        sensitivityRows(sensitivityCursor, :) = {vcal(iv), 'gamma_multiplier', multiplier, 0, ...
            params.gama * multiplier, finalYield, finalYield - propagate(detachment, baseCapture)};
    end
    for shift = energyShifts
        capture = demkovS14(vcal(iv), baseEnergy + shift, params.gama);
        finalYield = propagate(detachment, capture);
        sensitivityCursor = sensitivityCursor + 1;
        sensitivityRows(sensitivityCursor, :) = {vcal(iv), 'energy_defect_shift_au', 1, shift, ...
            params.gama, finalYield, finalYield - propagate(detachment, baseCapture)};
    end
end
sensitivity = cell2table(sensitivityRows, 'VariableNames', { ...
    'incident_velocity_au','sensitivity_axis','gamma_multiplier','energy_defect_shift_au', ...
    'gamma_used','final_yield','change_from_si_s14_baseline'});
writetable(sensitivity, fullfile(tableRoot, 'demkov_local_sensitivity.csv'));

auditSummary = struct();
auditSummary.schema_version = 1;
auditSummary.octopus_launched = false;
auditSummary.classical_replay_max_event_height_error_bohr = max(replay.max_event_height_error_bohr);
auditSummary.classical_replay_max_legacy_turning_height_error_bohr = max(replay.legacy_turning_height_absolute_error_bohr);
auditSummary.maximum_true_minus_saved_turning_height_bohr = max(abs(replay.true_minus_saved_turning_height_bohr));
auditSummary.si_s14_vs_archived_max_capture_probability_difference = baselineCaptureMaxError;
auditSummary.total_reconstructed_events = height(events);
auditSummary.td_database_covered_events = sum(events.within_td_database);
auditSummary.td_database_out_of_domain_events = sum(~events.within_td_database);
auditSummary.td_database_coverage_fraction = mean(events.within_td_database);
auditSummary.q20_points_aligned = height(q20Td);
auditSummary.experiment_points_in_PdE5 = size(readmatrix(fullfile(processedRoot, 'expe.csv')), 1);
auditSummary.formal_Wgeom_status = 'diagnostic_only: six deterministic trajectories, no experimental trajectory weights or impact-parameter ensemble';
auditSummary.demkov_provenance_status = 'SI S14 implemented exactly; archived zPcapture mismatch and gamma literature-to-number provenance are disclosed';
writeJson(fullfile(tableRoot, 'p0_nonproduction_analysis_summary.json'), auditSummary);

summary = auditSummary;
fprintf('P0 non-production analysis complete: %d events, %d Q20 points.\n', height(events), height(q20Td));
end

function path = findOne(root, name, requiredFragment)
if nargin < 3
    requiredFragment = '';
end
matches = dir(fullfile(root, '**', name));
if ~isempty(requiredFragment)
    matches = matches(contains({matches.folder}, requiredFragment));
end
if numel(matches) ~= 1
    error('Expected one match for %s under %s; found %d.', name, root, numel(matches));
end
path = fullfile(matches(1).folder, matches(1).name);
end

function values = readOctopusNumeric(path)
fid = fopen(path, 'r');
if fid < 0
    error('Cannot open Octopus series: %s', path);
end
cleanup = onCleanup(@() fclose(fid));
columnCount = 0;
while ~feof(fid)
    line = strtrim(fgetl(fid));
    if isempty(line) || startsWith(line, '#')
        continue;
    end
    columnCount = numel(sscanf(line, '%f'));
    break;
end
if columnCount == 0
    error('No numeric rows found in %s.', path);
end
frewind(fid);
format = repmat('%f', 1, columnCount);
parsed = textscan(fid, format, 'CommentStyle', '#', 'Delimiter', {' ', '\t'}, ...
    'MultipleDelimsAsOne', true, 'CollectOutput', true);
values = parsed{1};
end

function result = interpDescending(x, y, target)
[sortedX, order] = sort(x(:));
sortedY = y(order, :);
result = interp1(sortedX, sortedY, target, 'linear');
end

function value = pearson(x, y)
c = corrcoef(x, y);
value = c(1,2);
end

function capture = demkovS14(velocity, energyDefect, gamma)
% Exact Supporting Information Eq. (S14); no additional height envelope.
capture = 0.5 .* sech(pi .* (energyDefect + velocity^2/2) ./ (2*gamma*velocity)).^2;
if any(~isfinite(capture), 'all') || any(capture < 0 | capture > 0.5, 'all')
    error('SI Eq. (S14) produced an invalid probability; clipping is forbidden.');
end
end

function finalYield = propagate(detachment, capture)
probability = 0;
for index = 1:numel(detachment)
    probability = (1 - probability) * capture(index) + probability * (1 - detachment(index));
end
finalYield = probability;
end

function writeJson(path, value)
text = jsonencode(value, PrettyPrint=true);
fid = fopen(path, 'w', 'n', 'UTF-8');
if fid < 0
    error('Cannot write JSON file: %s', path);
end
cleanup = onCleanup(@() fclose(fid));
fwrite(fid, text, 'char');
fwrite(fid, newline, 'char');
end
