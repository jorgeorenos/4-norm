% Initialize an SO(9) optimizer population from the lowest estimated
% 4-norm condition numbers among Haar rotations. The selected rotations are
% canonicalized by initialize_condition_number after their raw scores rank
% them; these numerical estimates do not certify global minima.
scriptDirectory = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptDirectory));
addpath(fullfile(projectRoot, 'src', '4-norm_handle'), '-begin');
addpath(fullfile(projectRoot, 'src', 'helpers'), '-begin');

rng(42, 'twister');

n = 9;
numSamples = 1000;
numRetained = 100;
batchSize = 100;
numRandomStarts = 75;
options = struct('NumRandomStarts', numRandomStarts, 'MaxIterations', 1000, ...
    'MaxWorkingMemoryMB', 128, 'ScreenIterations', 2, 'NumFinalists', 4);
normHandle = @norm4_columns;
conditionNumberFn = @(rotations) evaluate_condition_number_batch( ...
    rotations, normHandle, options);

timer = tic;
[filteredRotations, filteredConditionNumbers, filteredIndices] = ...
    initialize_condition_number(n, numSamples, batchSize, ...
    conditionNumberFn, numRetained);
elapsedSeconds = toc(timer);

fprintf('Selection and canonicalization time: %.3f seconds\n', elapsedSeconds);
fprintf('Mean: %.12f | median: %.12f | standard deviation: %.12f\n', ...
    mean(filteredConditionNumbers), median(filteredConditionNumbers), ...
    std(filteredConditionNumbers));

