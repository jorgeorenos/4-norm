% Show how deterministic Pattern Search moves through canonical SO(2).
% The background and the displayed path use the same canonical SO(2) domain.
% The black path contains only accepted Pattern Search iterates from a
% deliberately difficult sample; segments crossing canonical boundaries break.
scriptDirectory = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptDirectory));
addpath(fullfile(projectRoot, 'src', '4-norm_handle'), '-begin');
addpath(fullfile(projectRoot, 'src', 'helpers'), '-begin');
addpath(fullfile(projectRoot, 'src', 'optimizer'), '-begin');

K = 2;
numSamples = 10000;
batchSize = 1000;
numRandomStarts = 75;
angleLimitDegrees = 45;
cellWidthDegrees = 0.25;
seeds = struct('Haar', 42, 'Starts', 43);

% P = I preserves the objective of the SO(2) cells example: kappa_4(Q).
P = eye(K);
normOptions = struct('MaxIterations', 1000, ...
    'NormTolerance', 1e-10, 'StationarityTolerance', 1e-8, ...
    'FeasibilityTolerance', 1e-12, 'MaxWorkingMemoryMB', 128, ...
    'ScreenIterations', 2, 'NumFinalists', 4);
normHandle = @norm4_columns;
banks = create_norm4_start_banks(K, numRandomStarts, seeds.Starts);

rng(seeds.Haar, 'twister');
canonicalAngles = zeros(numSamples, 1);
conditionNumbers = zeros(numSamples, 1);
startRotation = zeros(K);
startScreeningValue = -Inf;

timer = tic;
for batchStart = 1:batchSize:numSamples
    batchEnd = min(batchStart + batchSize - 1, numSamples);
    batchIndices = batchStart:batchEnd;
    count = numel(batchIndices);
    rotations = zeros(K, K, count);
    for page = 1:count
        rotations(:,:,page) = haar_so(K);
    end
    canonicalRotations = canonicalize_so2_signed_permutations(rotations);
    batchValues = compute_cond4_pq_fixed_starts( ...
        P, canonicalRotations, banks, normHandle, normOptions);
    conditionNumbers(batchIndices) = batchValues;
    canonicalAngles(batchIndices) = reshape(atan2(canonicalRotations(2,1,:), ...
        canonicalRotations(1,1,:)), [], 1);
    [batchMaximum, batchMaximumIndex] = max(batchValues);
    if batchMaximum > startScreeningValue
        startScreeningValue = batchMaximum;
        startRotation = canonicalRotations(:,:,batchMaximumIndex);
    end
end
mapSeconds = toc(timer);

assert(all(isfinite(conditionNumbers)) && all(conditionNumbers > 0));
assert(all(abs(canonicalAngles) <= pi/4 + 1e-12));
assert(norm(startRotation'*startRotation - eye(K), 'fro') < 1e-12);
assert(abs(det(startRotation) - 1) < 1e-12);

anglesDegrees = rad2deg(canonicalAngles);
cellMap = occupied_condition_cells(anglesDegrees, conditionNumbers, ...
    angleLimitDegrees, cellWidthDegrees);
[minimumCondition, minimumIndex] = min(conditionNumbers);

% Use the same fixed banks as the background so every trajectory point lies
% on the displayed finite-start numerical objective.
searchOptions = struct('MeshTolerance', 1e-4, ...
    'MeshExpansionFactor', 2, 'MeshContractionFactor', 0.5, ...
    'MaxIterations', 30, ...
    'MaxFunctionEvaluations', 200, 'Display', 'off', 'TrackHistory', true);
searchTimer = tic;
searchResult = pattern_search_cond4(P, startRotation, banks, normHandle, ...
    normOptions, searchOptions);
searchSeconds = toc(searchTimer);

historyRotations = reshape(searchResult.history.rotations, K, K, []);
historyCanonical = canonicalize_so2_signed_permutations(historyRotations);
historyAngles = rad2deg(reshape(atan2(historyCanonical(2,1,:), ...
    historyCanonical(1,1,:)), [], 1));
historyValues = searchResult.history.values(1,:).';
acceptedColumns = [1; find(searchResult.history.accepted(1,2:end)).' + 1];
acceptedAngles = historyAngles(acceptedColumns);
acceptedValues = historyValues(acceptedColumns);
finalRotation = canonicalize_so2_signed_permutations(searchResult.bestRotation);
finalAngle = rad2deg(atan2(finalRotation(2,1), finalRotation(1,1)));

referenceAngles = linspace(-pi/4, pi/4, 721);
phi = (0:1439)*(2*pi/1440);
directions = [cos(phi); sin(phi)];
X = directions ./ norm4_columns(directions);
referenceValues = zeros(size(referenceAngles));
for angleIndex = 1:numel(referenceAngles)
    angle = referenceAngles(angleIndex);
    Q = [cos(angle), -sin(angle); sin(angle), cos(angle)];
    referenceValues(angleIndex) = max(norm4_columns(Q*X))*max(norm4_columns(Q'*X));
end

figureHandle = figure('Name', 'Trayectoria Pattern Search en SO(2)', 'Color', 'w');
layout = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
colormap(figureHandle, turbo(256));
axesMap = nexttile;
plot(axesMap, rad2deg(referenceAngles), referenceValues, '-', ...
    'Color', [0.55, 0.55, 0.55], 'LineWidth', 1.1, ...
    'DisplayName', 'Referencia angular');
hold(axesMap, 'on');
scatter(axesMap, cellMap.angle, cellMap.values, 20, log10(cellMap.values), ...
    'filled', 'MarkerEdgeColor', 'none', 'DisplayName', 'Minimo por celda');
plot_canonical_so2_segments(axesMap, acceptedAngles, acceptedValues, ...
    angleLimitDegrees, [0, 0, 0], 1.4);
plot(axesMap, acceptedAngles, acceptedValues, 'ko', ...
    'MarkerFaceColor', 'w', 'MarkerSize', 5, ...
    'DisplayName', 'Iterados aceptados');
plot(axesMap, acceptedAngles(1), acceptedValues(1), 'ks', ...
    'MarkerFaceColor', 'y', 'MarkerSize', 9, 'DisplayName', 'Inicio');
plot(axesMap, finalAngle, searchResult.bestValue, 'rp', ...
    'MarkerFaceColor', 'r', 'MarkerSize', 11, 'DisplayName', 'Solucion final');
hold(axesMap, 'off');
colorLimits = log10([min(cellMap.values), max(cellMap.values)]);
clim(axesMap, colorLimits);
colorbarHandle = colorbar(axesMap);
ylabel(colorbarHandle, 'log_{10}(minimo estimado de \kappa_4(Q))');
grid(axesMap, 'on');
xlim(axesMap, [-angleLimitDegrees, angleLimitDegrees]);
ylim(axesMap, [0.99, sqrt(2) + 0.02]);
xlabel(axesMap, '\theta canonico (grados)');
ylabel(axesMap, '\kappa_4(Q)');
title(axesMap, {'Trayectoria aceptada', 'sobre el mapa de celdas'});
legend(axesMap, 'Location', 'northwest');

axesProgress = nexttile;
outerIterations = 0:searchResult.history.outerIterations;
plot(axesProgress, outerIterations, historyValues, 'o-', ...
    'Color', [0.10, 0.25, 0.70], 'LineWidth', 1.3, ...
    'MarkerFaceColor', [0.10, 0.25, 0.70]);
ylabel(axesProgress, '\kappa_4(Q) estimada');
grid(axesProgress, 'on');
xlabel(axesProgress, 'Iteracion exterior');
title(axesProgress, {'Objetivo estimado', 'durante Pattern Search'});

sgtitle(layout, sprintf(['SO(2) | %s Q canonicas | inicio %.6f | ', ...
    'final %.6f'], thousands_label(numSamples), ...
    searchResult.initialValues(1), searchResult.bestValue));

fprintf('\nTrayectoria Pattern Search SO(2)\n');
fprintf('Mapa de %d Q en %.3f segundos; Pattern Search en %.3f segundos.\n', ...
    numSamples, mapSeconds, searchSeconds);
fprintf('Inicio: %.12f | final: %.12f | iteraciones: %d | evaluaciones: %d.\n', ...
    searchResult.initialValues(1), searchResult.bestValue, ...
    searchResult.iterations(1), searchResult.functionCounts(1));
fprintf('Angulo final: %.6f grados | minimo muestral del mapa: %.12f.\n', ...
    finalAngle, minimumCondition);

function cellMap = occupied_condition_cells(angles, values, limitDegrees, cellWidthDegrees)
    edges = -limitDegrees:cellWidthDegrees:limitDegrees;
    centers = edges(1:end-1) + cellWidthDegrees/2;
    angleBin = discretize(angles, edges);
    valid = ~isnan(angleBin);
    minima = accumarray(angleBin(valid), values(valid), ...
        [numel(centers), 1], @min, NaN);
    occupied = find(isfinite(minima));
    cellMap = struct('angle', centers(occupied).', 'values', minima(occupied));
end

function plot_canonical_so2_segments(axesHandle, angles, values, limitDegrees, color, lineWidth)
    for pointIndex = 2:numel(angles)
        if abs(angles(pointIndex) - angles(pointIndex - 1)) <= limitDegrees
            plot(axesHandle, angles(pointIndex - 1:pointIndex), ...
                values(pointIndex - 1:pointIndex), '-', 'Color', color, ...
                'LineWidth', lineWidth, 'HandleVisibility', 'off');
        end
    end
end

function label = thousands_label(value)
    label = sprintf('%d', value);
    positions = numel(label)-2:-3:2;
    for position = positions
        label = [label(1:position-1), ',', label(position:end)]; %#ok<AGROW>
    end
end