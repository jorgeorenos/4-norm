% End-to-end SO(9) Pattern Search example for minimizing kappa_4(P*Q).
% Fixed forward/inverse start banks are used first to select 100 rotations
% from 10000 Haar draws and then to define a deterministic search objective.
% A progressive schedule keeps 25 trajectories after iteration 10 and 5
% after iteration 30, preserving each survivor's complete search state.
% The returned estimate is local and does not certify a global minimum.
scriptDirectory = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptDirectory));
addpath(fullfile(projectRoot, 'src', '4-norm_handle'), '-begin');
addpath(fullfile(projectRoot, 'src', 'helpers'), '-begin');
addpath(fullfile(projectRoot, 'src', 'optimizer'), '-begin');


K = 9;
numSamples = 10000;
numRetained = 100;
batchSize = 1500;

seeds = struct('Impact', 41, 'Haar', 42, ...
    'Screening', 43, 'Optimization', 44);

% In a fitted SVAR, SigmaE and P come from the estimated reduced-form
% covariance. This synthetic example gives P its own reproducible seed.
rng(seeds.Impact, 'twister');
gaussianMatrix = randn(K);
SigmaE = gaussianMatrix*gaussianMatrix' + 0.25*eye(K);
P = chol((SigmaE + SigmaE')/2, 'lower');

normHandle = @norm4_columns;
numScreeningStarts = 75;
screeningBanks = create_norm4_start_banks( ...
    K, numScreeningStarts, seeds.Screening);
screeningNormOptions = struct('MaxIterations', 1000, ...
    'NormTolerance', 1e-10, 'StationarityTolerance', 1e-8, ...
    'FeasibilityTolerance', 1e-12, 'MaxWorkingMemoryMB', 128, ...
    'ScreenIterations', 2, 'NumFinalists', 4);
conditionNumberFn = @(rotations) compute_cond4_pq_fixed_starts( ...
    P, rotations, screeningBanks, normHandle, screeningNormOptions);

% The fixed screening banks consume no RNG. Therefore this seed controls only
% the 10000 Haar rotations, independently of the norm-estimation vectors.
rng(seeds.Haar, 'twister');
selectionTimer = tic;
[initialRotations, screeningValues, selectedIndices] = ...
    initialize_condition_number(K, numSamples, batchSize, ...
    conditionNumberFn, numRetained);
selectionSeconds = toc(selectionTimer);

% Pattern Search receives a second pair of fixed banks. Every initial point
% and every Givens poll candidate is compared with these same vectors.
numOptimizationStarts = 75;
optimizationBanks = create_norm4_start_banks( ...
    K, numOptimizationStarts, seeds.Optimization);
optimizationNormOptions = struct('MaxIterations', 500, ...
    'NormTolerance', 1e-10, 'StationarityTolerance', 1e-8, ...
    'FeasibilityTolerance', 1e-12, 'MaxWorkingMemoryMB', 128, ...
    'ScreenIterations', 2, 'NumFinalists', 4);
searchOptions = struct('InitialMeshSize', 0.04, ...
    'MeshTolerance', 0.001, 'MaxMeshSize', 0.16, ...
    'MeshExpansionFactor', 2, 'MeshContractionFactor', 0.5, ...
    'FunctionTolerance', 1e-6, 'MaxIterations', 50, ...
    'MaxFunctionEvaluations', 3000, 'Display', 'off', ...
    'StageIterationLimits', [10, 30, 50], ...
    'StageRetainedCounts', [numRetained, 25, 5]);

searchTimer = tic;
searchResult = pattern_search_cond4(P, initialRotations, ...
    optimizationBanks, normHandle, optimizationNormOptions, searchOptions);
searchSeconds = toc(searchTimer);

QStar = searchResult.bestRotation;
kappa4Star = searchResult.bestValue;
impactMatrix = P*QStar;
initialBest = min(searchResult.initialValues);
relativeImprovement = (initialBest - kappa4Star)/initialBest;
orthogonalityError = norm(QStar'*QStar - eye(K), 'fro');
determinantError = abs(det(QStar) - 1);
covarianceError = norm(impactMatrix*impactMatrix' - SigmaE, 'fro') / ...
    norm(SigmaE, 'fro');

assert(isequal(size(initialRotations), [K, K, numRetained]));
assert(all(searchResult.finalValues <= searchResult.initialValues + 1e-12));
assert(orthogonalityError < 1e-12*max(1,K));
assert(determinantError < 1e-12*max(1,K));
assert(covarianceError < 1e-12*max(1,K));

fprintf('\nSO(9) Pattern Search example\n');
fprintf('Selection: %d of %d rotations in %.3f seconds.\n', ...
    numRetained, numSamples, selectionSeconds);
fprintf('Screening raw-estimate range: [%.12f, %.12f].\n', ...
    screeningValues(1), screeningValues(end));
fprintf('Pattern Search time: %.3f seconds.\n', searchSeconds);
fprintf(['Tracked objective evaluations across trajectories ', ...
    '(excluding final canonical reevaluation): %d.\n'], ...
    sum(searchResult.functionCounts));
fprintf('Progressive pruning: %d trajectories; full budget: %d.\n', ...
    sum(strcmp(searchResult.terminationReasons, 'StagePruned')), ...
    sum(strcmp(searchResult.terminationReasons, 'MaxIterations')));
fprintf('Best canonical initial estimate: %.12f.\n', initialBest);
fprintf('Best final estimate: %.12f.\n', kappa4Star);
fprintf('Relative estimated improvement: %.6f%%.\n', 100*relativeImprovement);
fprintf('Winning retained start: %d (Haar draw %d).\n', ...
    searchResult.bestStartIndex, selectedIndices(searchResult.bestStartIndex));
fprintf('Orthogonality error: %.3e | determinant error: %.3e.\n', ...
    orthogonalityError, determinantError);
fprintf('Relative covariance reconstruction error: %.3e.\n', covarianceError);
fprintf('Final canonical Q:\n');
disp(QStar);
