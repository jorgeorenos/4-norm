% Plot the winning full-budget Pattern Search path for a synthetic SO(9) case.
scriptDirectory = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptDirectory));
addpath(fullfile(projectRoot, 'src', '4-norm_handle'), '-begin');
addpath(fullfile(projectRoot, 'src', 'helpers'), '-begin');
addpath(fullfile(projectRoot, 'src', 'optimizer'), '-begin');

K = 9;
numSamples = 10000;
numPatternSearchStarts = 100;
batchSize = 1500;
numNormStarts = 75;
seeds = struct('Impact', 41, 'Haar', 42, 'Starts', 43);

rng(seeds.Impact, 'twister');
gaussianMatrix = randn(K);
SigmaE = gaussianMatrix*gaussianMatrix' + 0.25*eye(K);
P = chol((SigmaE + SigmaE')/2, 'lower');

normHandle = @norm4_columns;
normOptions = struct('MaxIterations', 500, ...
    'NormTolerance', 1e-10, 'StationarityTolerance', 1e-8, ...
    'FeasibilityTolerance', 1e-12, 'MaxWorkingMemoryMB', 128, ...
    'ScreenIterations', 2, 'NumFinalists', 4);
banks = create_norm4_start_banks(K, numNormStarts, seeds.Starts);
conditionNumberFn = @(rotations) compute_cond4_pq_fixed_starts( ...
    P, canonicalize_rotations_by_trace_batch(rotations), ...
    banks, normHandle, normOptions);

rng(seeds.Haar, 'twister');
selectionTimer = tic;
[initialRotations, screeningValues, selectedIndices] = ...
    initialize_condition_number(K, numSamples, batchSize, ...
    conditionNumberFn, numPatternSearchStarts);
selectionSeconds = toc(selectionTimer);
assert(size(initialRotations, 3) == numPatternSearchStarts);
assert(issorted(screeningValues));

baselineValue = compute_cond4_pq_fixed_starts( ...
    P, eye(K), banks, normHandle, normOptions);
searchOptions = struct('MeshTolerance', 0.001, ...
    'MeshExpansionFactor', 2, 'MeshContractionFactor', 0.5, ...
    'MaxIterations', 50, 'Display', 'off', 'TrackHistory', true);
searchTimer = tic;
searchResult = pattern_search_cond4(P, initialRotations, banks, ...
    normHandle, normOptions, searchOptions);
searchSeconds = toc(searchTimer);
meshStoppedCount = sum(strcmp(searchResult.terminationReasons, 'MeshTolerance'));
roundLimitCount = sum(strcmp(searchResult.terminationReasons, 'MaxIterations'));
assert(meshStoppedCount + roundLimitCount == numPatternSearchStarts);

bestStartIndex = searchResult.bestStartIndex;
bestTrajectoryIterations = searchResult.iterations(bestStartIndex);
historyValues = searchResult.history.values( ...
    bestStartIndex, 1:bestTrajectoryIterations + 1).';
outerIterations = (0:bestTrajectoryIterations).';
assert(numel(searchResult.initialValues) == numPatternSearchStarts);
assert(abs(historyValues(1) - searchResult.initialValues(bestStartIndex)) < 1e-12);
assert(abs(historyValues(end) - searchResult.bestValue) < 1e-12);
assert(all(diff(historyValues) <= 1e-12));
assert(all(searchResult.finalValues <= searchResult.initialValues + 1e-12));

bestRotation = searchResult.bestRotation;
impactMatrix = P*bestRotation;
orthogonalityError = norm(bestRotation'*bestRotation - eye(K), 'fro');
determinantError = abs(det(bestRotation) - 1);
covarianceError = norm(impactMatrix*impactMatrix' - SigmaE, 'fro') / ...
    norm(SigmaE, 'fro');
assert(orthogonalityError < 1e-12*max(1,K));
assert(determinantError < 1e-12*max(1,K));
assert(covarianceError < 1e-12*max(1,K));

figureHandle = figure('Name', 'Pattern Search ganador en SO(9)', ...
    'Color', 'w', 'Position', [100, 100, 1000, 600]);
axesHandle = axes('Parent', figureHandle);
plot(axesHandle, outerIterations, historyValues, 'o-', ...
    'Color', [0.10, 0.25, 0.70], 'LineWidth', 1.4, ...
    'MarkerFaceColor', [0.10, 0.25, 0.70], ...
    'DisplayName', 'Trayectoria ganadora entre 100 inicios');
hold(axesHandle, 'on');
yline(axesHandle, baselineValue, 'k:', 'LineWidth', 1.5, ...
    'DisplayName', 'Referencia: kappa_4(P), Q = I');
hold(axesHandle, 'off');
grid(axesHandle, 'on');
xlabel(axesHandle, 'iteracion');
ylabel(axesHandle, 'kappa_4(PQ) estimada');
title(axesHandle, sprintf(['SO(9): ganadora %d (muestra %d), ', ...
    '%.6f a %.6f'], bestStartIndex, selectedIndices(bestStartIndex), ...
    searchResult.initialValues(bestStartIndex), searchResult.bestValue));
progressYLimits = ylim(axesHandle);
axisPadding = 0.05*diff(progressYLimits);
ylim(axesHandle, [progressYLimits(1) - axisPadding, ...
    progressYLimits(2) + axisPadding]);
legend(axesHandle, 'Location', 'best');

outputDirectory = fullfile(projectRoot, 'docs', 'patter_search');
if ~isfolder(outputDirectory)
    mkdir(outputDirectory);
end
exportgraphics(figureHandle, ...
    fullfile(outputDirectory, 'pattern_search_cond4_9x9.png'), ...
    'Resolution', 200);

fprintf('\nPattern Search SO(9)\n');
fprintf('Selection: %d lowest estimates from %d Haar rotations in %.3f seconds.\n', ...
    numPatternSearchStarts, numSamples, selectionSeconds);
fprintf('Pattern Search: %.3f seconds; all 100 trajectories used up to 50 rounds.\n', ...
    searchSeconds);
fprintf('Stops: %d at mesh tolerance; %d at the round limit; no pruning.\n', ...
    meshStoppedCount, roundLimitCount);
fprintf(['Winner %d (Haar draw %d): start %.12f | final %.12f | ', ...
    'rounds %d | evaluations %d.\n'], bestStartIndex, ...
    selectedIndices(bestStartIndex), searchResult.initialValues(bestStartIndex), ...
    searchResult.bestValue, bestTrajectoryIterations, ...
    searchResult.functionCounts(bestStartIndex));
fprintf('Total objective evaluations across 100 trajectories: %d.\n', ...
    sum(searchResult.functionCounts));
fprintf('Reference kappa_4(P), Q=I: %.12f.\n', baselineValue);
fprintf('Orthogonality %.3e | determinant %.3e | covariance %.3e.\n', ...
    orthogonalityError, determinantError, covarianceError);
