% Plot the winning full-budget Pattern Search path for a synthetic SO(3) case.
scriptDirectory = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptDirectory));
addpath(fullfile(projectRoot, 'src', '4-norm_handle'), '-begin');
addpath(fullfile(projectRoot, 'src', 'helpers'), '-begin');
addpath(fullfile(projectRoot, 'src', 'optimizer'), '-begin');

K = 3;
numSamples = 10000;
numPatternSearchStarts = 100;
batchSize = 100;
numNormStarts = 75;
seeds = struct('Impact', 41, 'Haar', 42, 'Starts', 43);

rng(seeds.Impact, 'twister');
gaussianMatrix = randn(K);
SigmaE = gaussianMatrix*gaussianMatrix' + 0.25*eye(K);
P = chol((SigmaE + SigmaE')/2, 'lower');
normOptions = struct('MaxIterations', 1000, ...
    'NormTolerance', 1e-10, 'StationarityTolerance', 1e-8, ...
    'FeasibilityTolerance', 1e-12, 'MaxWorkingMemoryMB', 128, ...
    'ScreenIterations', 2, 'NumFinalists', 4);
normHandle = @norm4_columns;
banks = create_norm4_start_banks(K, numNormStarts, seeds.Starts);

rng(seeds.Haar, 'twister');
conditionNumbers = zeros(numSamples, 1);
sampleRotations = zeros(K, K, numSamples);

selectionTimer = tic;
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
    sampleRotations(:,:,batchIndices) = canonicalRotations;
end
selectionSeconds = toc(selectionTimer);

assert(all(isfinite(conditionNumbers)) && all(conditionNumbers > 0));
[~, sampleOrder] = sortrows([conditionNumbers, (1:numSamples).'], [1, 2]);
selectedIndices = sampleOrder(1:numPatternSearchStarts);
initialRotations = sampleRotations(:,:,selectedIndices);
selectedSampleValues = conditionNumbers(selectedIndices);
assert(issorted(selectedSampleValues));
assert(size(initialRotations, 3) == numPatternSearchStarts);

searchOptions = struct('MeshTolerance', 0.001, ...
    'MeshExpansionFactor', 2, 'MeshContractionFactor', 0.5, ...
    'MaxIterations', 40, ...
    'Display', 'off', 'TrackHistory', true);
searchTimer = tic;
searchResult = pattern_search_cond4(P, initialRotations, banks, normHandle, ...
    normOptions, searchOptions);
searchSeconds = toc(searchTimer);
bestStartIndex = searchResult.bestStartIndex;
bestTrajectoryIterations = searchResult.iterations(bestStartIndex);
assert(numel(searchResult.initialValues) == numPatternSearchStarts);
baselineValue = compute_cond4_pq_fixed_starts( ...
    P, eye(K), banks, normHandle, normOptions);

historyValues = searchResult.history.values( ...
    bestStartIndex,1:bestTrajectoryIterations + 1).';

figureHandle = figure('Name', 'Pattern Search ganador en SO(3)', ...
    'Color', 'w', 'Position', [100, 100, 1000, 600]);
axesProgress = axes('Parent', figureHandle);
outerIterations = (0:bestTrajectoryIterations).';
plot(axesProgress, outerIterations, historyValues, 'o-', ...
    'Color', [0.10, 0.25, 0.70], 'LineWidth', 1.4, ...
    'MarkerFaceColor', [0.10, 0.25, 0.70], ...
    'DisplayName', 'Trayectoria ganadora entre 100 inicios');
hold(axesProgress, 'on');
yline(axesProgress, baselineValue, 'k:', 'LineWidth', 1.5, ...
    'DisplayName', 'Referencia: kappa_4(P), Q = I');
hold(axesProgress, 'off');
grid(axesProgress, 'on');
xlabel(axesProgress, 'iteracion');
ylabel(axesProgress, 'kappa_4(PQ) estimada');
title(axesProgress, sprintf(['SO(3): ganadora %d (muestra %d), ', ...
    '%.6f a %.6f'], bestStartIndex, selectedIndices(bestStartIndex), ...
    searchResult.initialValues(bestStartIndex), searchResult.bestValue));
progressYLimits = ylim(axesProgress);
axisPadding = 0.05*diff(progressYLimits);
ylim(axesProgress, [progressYLimits(1) - axisPadding, ...
    progressYLimits(2) + axisPadding]);
legend(axesProgress, 'Location', 'best');

outputDirectory = fullfile(projectRoot, 'docs', 'patter_search');
if ~isfolder(outputDirectory)
    mkdir(outputDirectory);
end
exportgraphics(figureHandle, ...
    fullfile(outputDirectory, 'pattern_search_cond4_3x3.png'), ...
    'Resolution', 200);

fprintf('\nTrayectoria Pattern Search SO(3)\n');
fprintf('Selección de %d Q en %.3f segundos; Pattern Search en %.3f segundos.\n', ...
    numSamples, selectionSeconds, searchSeconds);
fprintf('Se refinaron las %d Q de menor estimacion entre las %d muestreadas.\n', ...
    numPatternSearchStarts, numSamples);
fprintf(['Trayectoria ganadora (%d): inicio %.12f | final %.12f | ', ...
    'rondas %d | evaluaciones %d.\n'], bestStartIndex, ...
    searchResult.initialValues(bestStartIndex), searchResult.bestValue, ...
    bestTrajectoryIterations, searchResult.functionCounts(bestStartIndex));
fprintf('Evaluaciones totales en las 100 trayectorias: %d.\n', ...
    sum(searchResult.functionCounts));
fprintf('Referencia kappa_4(P), Q=I: %.12f.\n', baselineValue);