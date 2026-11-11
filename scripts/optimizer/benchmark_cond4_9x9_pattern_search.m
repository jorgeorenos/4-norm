% Reproducible benchmark for SO(9) Pattern Search configurations.
% The condition-number estimator and its fixed banks remain identical across
% configurations. Three repetitions run in alternating order; reported times
% are medians.
scriptDirectory = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptDirectory));
addpath(fullfile(projectRoot, 'src', '4-norm_handle'), '-begin');
addpath(fullfile(projectRoot, 'src', 'helpers'), '-begin');
addpath(fullfile(projectRoot, 'src', 'optimizer'), '-begin');

K = 9;
numSamples = 10000;
numRetained = 100;
batchSize = 1500;
numBenchmarkRepeats = 3;

seeds = struct('Impact', 41, 'Haar', 42, ...
    'Screening', 43, 'Optimization', 44);

rng(seeds.Impact, 'twister');
gaussianMatrix = randn(K);
SigmaE = gaussianMatrix*gaussianMatrix' + 0.25*eye(K);
P = chol((SigmaE + SigmaE')/2, 'lower');

normHandle = @norm4_columns;
screeningBanks = create_norm4_start_banks(K, 75, seeds.Screening);
screeningNormOptions = struct('MaxIterations', 1000, ...
    'NormTolerance', 1e-10, 'StationarityTolerance', 1e-8, ...
    'FeasibilityTolerance', 1e-12, 'MaxWorkingMemoryMB', 128, ...
    'ScreenIterations', 2, 'NumFinalists', 4);
conditionNumberFn = @(rotations) compute_cond4_pq_fixed_starts( ...
    P, rotations, screeningBanks, normHandle, screeningNormOptions);

rng(seeds.Haar, 'twister');
selectionTimer = tic;
[initialRotations, ~, selectedIndices] = initialize_condition_number( ...
    K, numSamples, batchSize, conditionNumberFn, numRetained);
selectionSeconds = toc(selectionTimer);

optimizationBanks = create_norm4_start_banks(K, 75, seeds.Optimization);
optimizationNormOptions = struct('MaxIterations', 500, ...
    'NormTolerance', 1e-10, 'StationarityTolerance', 1e-8, ...
    'FeasibilityTolerance', 1e-12, 'MaxWorkingMemoryMB', 128, ...
    'ScreenIterations', 2, 'NumFinalists', 4);
baseSearchOptions = struct('MeshTolerance', 0.001, ...
    'MeshExpansionFactor', 2, 'MeshContractionFactor', 0.5, ...
    'MaxIterations', 50, ...
    'Display', 'off');

configurationNames = {'Baseline 100', 'Configured 50', 'Progressive 100-25-5'};
configurationStartCounts = [100, 50, 100];
configurationStageLimits = {[], [], [10, 30, 50]};
configurationRetainedCounts = {[], [], [100, 25, 5]};
configurationCount = numel(configurationNames);

% Warm up the execution path without including it in a measured configuration.
warmupOptions = baseSearchOptions;
warmupOptions.MaxIterations = 1;
pattern_search_cond4(P, initialRotations(:,:,1), optimizationBanks, ...
    normHandle, optimizationNormOptions, warmupOptions);

measuredSeconds = zeros(configurationCount, numBenchmarkRepeats);
configurationResults = cell(configurationCount, 1);
for repetition = 1:numBenchmarkRepeats
    if mod(repetition, 2) == 1
        executionOrder = 1:configurationCount;
    else
        executionOrder = configurationCount:-1:1;
    end
    for orderIndex = 1:configurationCount
        configurationIndex = executionOrder(orderIndex);
        startCount = configurationStartCounts(configurationIndex);
        searchOptions = baseSearchOptions;
        if ~isempty(configurationStageLimits{configurationIndex})
            searchOptions.StageIterationLimits = ...
                configurationStageLimits{configurationIndex};
            searchOptions.StageRetainedCounts = ...
                configurationRetainedCounts{configurationIndex};
        end
        searchTimer = tic;
        searchResult = pattern_search_cond4(P, ...
            initialRotations(:,:,1:startCount), optimizationBanks, ...
            normHandle, optimizationNormOptions, searchOptions);
        measuredSeconds(configurationIndex, repetition) = toc(searchTimer);
        if isempty(configurationResults{configurationIndex})
            configurationResults{configurationIndex} = searchResult;
        else
            assert(isequal(searchResult.bestValue, ...
                configurationResults{configurationIndex}.bestValue));
            assert(isequal(searchResult.functionCounts, ...
                configurationResults{configurationIndex}.functionCounts));
        end
    end
end

benchmarkResults = repmat(struct('name', '', 'searchSeconds', 0, ...
    'totalSeconds', 0, 'bestValue', 0, 'bestStartIndex', 0, ...
    'haarDrawIndex', 0, 'functionCount', 0, 'stagePrunedCount', 0, ...
    'fullBudgetCount', 0, 'timeReductionPercent', 0, ...
    'valueDifferencePercent', 0), configurationCount, 1);
for configurationIndex = 1:configurationCount
    searchResult = configurationResults{configurationIndex};
    QCandidate = searchResult.bestRotation;
    impactMatrix = P*QCandidate;
    assert(norm(QCandidate'*QCandidate - eye(K), 'fro') < 1e-12*max(1,K));
    assert(abs(det(QCandidate) - 1) < 1e-12*max(1,K));
    assert(norm(impactMatrix*impactMatrix' - SigmaE, 'fro') / ...
        norm(SigmaE, 'fro') < 1e-12*max(1,K));

    benchmarkResults(configurationIndex).name = ...
        configurationNames{configurationIndex};
    benchmarkResults(configurationIndex).searchSeconds = ...
        median(measuredSeconds(configurationIndex,:));
    benchmarkResults(configurationIndex).totalSeconds = selectionSeconds + ...
        benchmarkResults(configurationIndex).searchSeconds;
    benchmarkResults(configurationIndex).bestValue = searchResult.bestValue;
    benchmarkResults(configurationIndex).bestStartIndex = ...
        searchResult.bestStartIndex;
    benchmarkResults(configurationIndex).haarDrawIndex = ...
        selectedIndices(searchResult.bestStartIndex);
    benchmarkResults(configurationIndex).functionCount = ...
        sum(searchResult.functionCounts);
    benchmarkResults(configurationIndex).stagePrunedCount = ...
        sum(strcmp(searchResult.terminationReasons, 'StagePruned'));
    benchmarkResults(configurationIndex).fullBudgetCount = ...
        sum(strcmp(searchResult.terminationReasons, 'MaxIterations'));
end

baselineSeconds = benchmarkResults(1).searchSeconds;
baselineValue = benchmarkResults(1).bestValue;
for configurationIndex = 1:configurationCount
    benchmarkResults(configurationIndex).timeReductionPercent = 100 * ...
        (baselineSeconds - benchmarkResults(configurationIndex).searchSeconds) / ...
        baselineSeconds;
    benchmarkResults(configurationIndex).valueDifferencePercent = 100 * ...
        (benchmarkResults(configurationIndex).bestValue - baselineValue) / ...
        baselineValue;
end

fprintf('\nSO(9) Pattern Search reproducible benchmark\n');
fprintf('Selection time shared by all configurations: %.3f seconds.\n', ...
    selectionSeconds);
fprintf('%-24s %10s %10s %16s %12s %12s\n', ...
    'Configuration', 'Search(s)', 'Reduction', 'Best kappa4', ...
    'Tracked evals', 'Value gap');
for configurationIndex = 1:configurationCount
    entry = benchmarkResults(configurationIndex);
    fprintf('%-24s %10.3f %9.2f%% %16.12f %12d %11.4f%%\n', ...
        entry.name, entry.searchSeconds, entry.timeReductionPercent, ...
        entry.bestValue, entry.functionCount, entry.valueDifferencePercent);
end
