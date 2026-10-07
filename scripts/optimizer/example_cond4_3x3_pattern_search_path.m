% Show how deterministic Pattern Search moves through canonical SO(3).
% The background and displayed path use the same canonical SO(3) domain.
% The path contains only accepted Pattern Search iterates from a high-value
% sampled rotation; visual segments break at representative discontinuities.
scriptDirectory = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptDirectory));
addpath(fullfile(projectRoot, 'src', '4-norm_handle'), '-begin');
addpath(fullfile(projectRoot, 'src', 'helpers'), '-begin');
addpath(fullfile(projectRoot, 'src', 'optimizer'), '-begin');

K = 3;
numSamples = 10000;
batchSize = 100;
numRandomStarts = 75;
angleLimitDegrees = 50;
cellWidthDegrees = 1.755;
maximumDisplayedCells = 15000;
seeds = struct('Impact', 41, 'Haar', 42, 'Starts', 43, 'Display', 44);

rng(seeds.Impact, 'twister');
gaussianMatrix = randn(K);
SigmaE = gaussianMatrix*gaussianMatrix' + 0.25*eye(K);
P = chol((SigmaE + SigmaE')/2, 'lower');
normOptions = struct('MaxIterations', 1000, ...
    'NormTolerance', 1e-10, 'StationarityTolerance', 1e-8, ...
    'FeasibilityTolerance', 1e-12, 'MaxWorkingMemoryMB', 128, ...
    'ScreenIterations', 2, 'NumFinalists', 4);
normHandle = @norm4_columns;
banks = create_norm4_start_banks(K, numRandomStarts, seeds.Starts);

rng(seeds.Haar, 'twister');
angles = zeros(3, numSamples);
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
    canonicalRotations = canonicalize_so3_signed_permutations(rotations);
    batchValues = compute_cond4_pq_fixed_starts( ...
        P, canonicalRotations, banks, normHandle, normOptions);
    conditionNumbers(batchIndices) = batchValues;
    angles(:,batchIndices) = so3_zyx_euler_angles(canonicalRotations);
    [batchMaximum, batchMaximumIndex] = max(batchValues);
    if batchMaximum > startScreeningValue
        startScreeningValue = batchMaximum;
        startRotation = canonicalRotations(:,:,batchMaximumIndex);
    end
end
mapSeconds = toc(timer);

assert(all(isfinite(conditionNumbers)) && all(conditionNumbers > 0));
assert(all(isfinite(angles(:))));
assert(norm(startRotation'*startRotation - eye(K), 'fro') < 1e-12*max(1,K));
assert(abs(det(startRotation) - 1) < 1e-12*max(1,K));

anglesDegrees = rad2deg(angles);
alpha = anglesDegrees(1,:).';
beta = anglesDegrees(2,:).';
theta = anglesDegrees(3,:).';
cellMap = occupied_condition_cells(alpha, beta, theta, conditionNumbers, ...
    angleLimitDegrees, cellWidthDegrees);
occupiedCount = numel(cellMap.values);
displayedCount = min(maximumDisplayedCells, occupiedCount);
displayStream = RandStream('mt19937ar', 'Seed', seeds.Display);
displayPositions = randperm(displayStream, occupiedCount, displayedCount);

% The map and Pattern Search share the same fixed banks so the path is on
% precisely the numerical objective encoded by the cell colors.
searchOptions = struct('InitialMeshSize', 0.16, ...
    'MeshTolerance', 0.001, 'MaxMeshSize', 0.32, ...
    'MeshExpansionFactor', 2, 'MeshContractionFactor', 0.5, ...
    'FunctionTolerance', 1e-6, 'MaxIterations', 40, ...
    'MaxFunctionEvaluations', 400, 'Display', 'off', 'TrackHistory', true);
searchTimer = tic;
searchResult = pattern_search_cond4(P, startRotation, banks, normHandle, ...
    normOptions, searchOptions);
searchSeconds = toc(searchTimer);

historyRotations = reshape(searchResult.history.rotations, K, K, []);
historyCanonical = canonicalize_so3_signed_permutations(historyRotations);
historyAngles = rad2deg(so3_zyx_euler_angles(historyCanonical));
historyValues = searchResult.history.values(1,:).';
acceptedColumns = [1; find(searchResult.history.accepted(1,2:end)).' + 1];
acceptedAngles = historyAngles(:,acceptedColumns);
acceptedValues = historyValues(acceptedColumns);
finalRotation = canonicalize_so3_signed_permutations(searchResult.bestRotation);
finalAngles = rad2deg(so3_zyx_euler_angles(finalRotation));

figureHandle = figure('Name', 'Trayectoria Pattern Search en SO(3)', 'Color', 'w');
layout = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
colormap(figureHandle, turbo(256));
axesMap = nexttile;
scatter3(axesMap, cellMap.alpha(displayPositions), ...
    cellMap.beta(displayPositions), cellMap.theta(displayPositions), 10, ...
    log10(cellMap.values(displayPositions)), 'filled', 'MarkerEdgeColor', 'none', ...
    'MarkerFaceAlpha', 0.45);
hold(axesMap, 'on');
plot_canonical_so3_segments(axesMap, acceptedAngles, 30, [1, 1, 1], 5);
plot_canonical_so3_segments(axesMap, acceptedAngles, 30, [0, 0, 0], 2);
plot3(axesMap, acceptedAngles(1,:), acceptedAngles(2,:), acceptedAngles(3,:), ...
    'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'LineWidth', 1.2, ...
    'DisplayName', 'Iterados aceptados');
plot3(axesMap, acceptedAngles(1,1), acceptedAngles(2,1), acceptedAngles(3,1), ...
    'ks', 'MarkerFaceColor', 'y', 'MarkerSize', 11, 'LineWidth', 1.2, ...
    'DisplayName', 'Inicio');
plot3(axesMap, finalAngles(1), finalAngles(2), finalAngles(3), 'rp', ...
    'MarkerFaceColor', 'r', 'MarkerSize', 15, 'LineWidth', 1.2, ...
    'DisplayName', 'Solucion final');
hold(axesMap, 'off');
colorLimits = log10([min(cellMap.values), max(cellMap.values)]);
clim(axesMap, colorLimits);
colorbarHandle = colorbar(axesMap);
ylabel(colorbarHandle, 'log_{10}(minimo estimado de \kappa_4(PQ))');
axis(axesMap, 'vis3d');
xlim(axesMap, [-angleLimitDegrees, angleLimitDegrees]);
ylim(axesMap, [-angleLimitDegrees, angleLimitDegrees]);
zlim(axesMap, [-angleLimitDegrees, angleLimitDegrees]);
grid(axesMap, 'on');
box(axesMap, 'on');
view(axesMap, [-38, 22]);
xlabel(axesMap, '\alpha (grados)');
ylabel(axesMap, '\beta (grados)');
zlabel(axesMap, '\theta (grados)');
title(axesMap, {sprintf('%s/%s celdas ocupadas', ...
    thousands_label(displayedCount), thousands_label(occupiedCount)), ...
    'Negro: iterados', 'Amarillo: inicio | rojo: final'});

axesProgress = nexttile;
outerIterations = 0:searchResult.history.outerIterations;
plot(axesProgress, outerIterations, historyValues, 'o-', ...
    'Color', [0.10, 0.25, 0.70], 'LineWidth', 1.3, ...
    'MarkerFaceColor', [0.10, 0.25, 0.70]);
ylabel(axesProgress, '\kappa_4(PQ) estimada');
grid(axesProgress, 'on');
xlabel(axesProgress, 'Iteracion exterior');
title(axesProgress, {'Objetivo estimado', 'durante Pattern Search'});

sgtitle(layout, sprintf(['SO(3) | %s Q canonicas | inicio %.6f | ', ...
    'final %.6f'], thousands_label(numSamples), ...
    searchResult.initialValues(1), searchResult.bestValue));

fprintf('\nTrayectoria Pattern Search SO(3)\n');
fprintf('Mapa de %d Q en %.3f segundos; Pattern Search en %.3f segundos.\n', ...
    numSamples, mapSeconds, searchSeconds);
fprintf('Inicio: %.12f | final: %.12f | iteraciones: %d | evaluaciones: %d.\n', ...
    searchResult.initialValues(1), searchResult.bestValue, ...
    searchResult.iterations(1), searchResult.functionCounts(1));
fprintf('Iterados aceptados: %d | valor de inicio seleccionado: %.12f.\n', ...
    numel(acceptedValues), startScreeningValue);

function cellMap = occupied_condition_cells(alpha, beta, theta, values, limitDegrees, cellWidthDegrees)
    edges = -limitDegrees:cellWidthDegrees:limitDegrees;
    centers = edges(1:end-1) + cellWidthDegrees/2;
    binCount = numel(centers);
    alphaBin = discretize(alpha, edges);
    betaBin = discretize(beta, edges);
    thetaBin = discretize(theta, edges);
    valid = ~isnan(alphaBin) & ~isnan(betaBin) & ~isnan(thetaBin);
    linearBins = sub2ind([binCount, binCount, binCount], alphaBin(valid), ...
        betaBin(valid), thetaBin(valid));
    minima = accumarray(linearBins, values(valid), [binCount^3, 1], @min, NaN);
    occupied = find(isfinite(minima));
    [alphaIndex, betaIndex, thetaIndex] = ind2sub( ...
        [binCount, binCount, binCount], occupied);
    cellMap = struct('alpha', centers(alphaIndex).', ...
        'beta', centers(betaIndex).', 'theta', centers(thetaIndex).', ...
        'values', minima(occupied));
end

function plot_canonical_so3_segments(axesHandle, angles, maximumStepDegrees, color, lineWidth)
    for pointIndex = 2:size(angles, 2)
        if norm(angles(:,pointIndex) - angles(:,pointIndex - 1)) <= maximumStepDegrees
            plot3(axesHandle, angles(1,pointIndex - 1:pointIndex), ...
                angles(2,pointIndex - 1:pointIndex), ...
                angles(3,pointIndex - 1:pointIndex), '-', 'Color', color, ...
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