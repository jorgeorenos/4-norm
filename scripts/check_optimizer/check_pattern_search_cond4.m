% Basic checks for deterministic Pattern Search of kappa_4(P*Q).
scriptDirectory = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptDirectory));
addpath(fullfile(projectRoot, 'src', '4-norm_handle'), '-begin');
addpath(fullfile(projectRoot, 'src', 'helpers'), '-begin');
addpath(fullfile(projectRoot, 'src', 'optimizer'), '-begin');

callerState = rng;
banks = create_norm4_start_banks(2, 32, 31415);
assert(isequal(rng, callerState));
assert(isequal(size(banks.forward), [2, 32]));
assert(isequal(size(banks.inverse), [2, 32]));
assert(all(isfinite(banks.forward(:))) && all(isfinite(banks.inverse(:))));
assert(all(any(banks.forward ~= 0, 1)) && all(any(banks.inverse ~= 0, 1)));
assert(isequal(banks, create_norm4_start_banks(2, 32, 31415)));

normOptions = struct('MaxIterations', 300, 'NormTolerance', 1e-12, ...
    'StationarityTolerance', 1e-10, 'FeasibilityTolerance', 1e-12, ...
    'MaxWorkingMemoryMB', 64, 'ScreenIterations', 0, 'NumFinalists', 4);
identityValue = compute_cond4_pq_fixed_starts( ...
    eye(2), eye(2), banks, @norm4_columns, normOptions);
assert(abs(identityValue - 1) < 1e-12);

P = diag([2, 1]);
diagonalValue = compute_cond4_pq_fixed_starts( ...
    P, eye(2), banks, @norm4_columns, normOptions);
assert(abs(diagonalValue - 2) < 1e-8);

Q45 = [1, -1; 1, 1] / sqrt(2);
Qbatch = cat(3, eye(2), Q45);
batchValues = compute_cond4_pq_fixed_starts( ...
    P, Qbatch, banks, @norm4_columns, normOptions);
scalarValues = [ ...
    compute_cond4_pq_fixed_starts(P, eye(2), banks, @norm4_columns, normOptions); ...
    compute_cond4_pq_fixed_starts(P, Q45, banks, @norm4_columns, normOptions)];
assert(max(abs(batchValues - scalarValues)) < 1e-14);

searchOptions = struct('InitialMeshSize', pi/8, 'MeshTolerance', 1e-4, ...
    'MaxMeshSize', pi/4, 'MeshExpansionFactor', 2, ...
    'MeshContractionFactor', 0.5, ...
    'MaxIterations', 40, 'Display', 'off', ...
    'TrackHistory', true);
result = pattern_search_cond4(eye(2), Q45, banks, @norm4_columns, ...
    normOptions, searchOptions);
assert(result.bestValue < result.initialValues(1) - 1e-4);
assert(abs(result.bestValue - 1) < 1e-6);
assert(norm(result.bestRotation'*result.bestRotation - eye(2), 'fro') < 1e-12);
assert(abs(det(result.bestRotation) - 1) < 1e-12);
assert(isequal(rng, callerState));
assert(result.functionCounts(1) >= 1);
assert(result.iterations(1) >= 1);
assert(isfield(result, 'terminationReasons'));
assert(result.history.outerIterations == result.iterations(1));
assert(size(result.history.rotations, 4) == result.iterations(1) + 1);
assert(result.history.values(1,1) == result.initialValues(1));
assert(result.history.values(1,end) == result.rawFinalValues(1));
assert(isequal(result.finalRotations, result.rawFinalRotations));
assert(isequal(result.finalValues, result.rawFinalValues));
assert(isequal(result.bestRotation, ...
    result.rawFinalRotations(:,:,result.bestStartIndex)));
assert(result.bestValue == result.rawFinalValues(result.bestStartIndex));
assert(~isfield(result.options, 'FunctionTolerance'));
assert(~isfield(result.options, 'MaxFunctionEvaluations'));
removedBudgetOptions = searchOptions;
removedBudgetOptions.MaxFunctionEvaluations = 1;
caughtRemovedBudget = false;
try
    pattern_search_cond4(eye(2), Q45, banks, @norm4_columns, ...
        normOptions, removedBudgetOptions);
catch exception
    caughtRemovedBudget = strcmp(exception.identifier, ...
        'optimizer:pattern_search_cond4:UnknownOption');
end
assert(caughtRemovedBudget);
acceptedSteps = find(result.history.accepted(1,2:end)) + 1;
assert(all(result.history.values(1,acceptedSteps) < ...
    result.history.values(1,acceptedSteps - 1)));
rejectedOptions = searchOptions;
rejectedOptions.FunctionTolerance = 0;
caughtRemovedOption = false;
try
    pattern_search_cond4(eye(2), Q45, banks, @norm4_columns, ...
        normOptions, rejectedOptions);
catch exception
    caughtRemovedOption = strcmp(exception.identifier, ...
        'optimizer:pattern_search_cond4:UnknownOption');
end
assert(caughtRemovedOption);

% Omitting mesh sizes uses the optimizer defaults, not script-specific sizes.
defaultMeshOptions = rmfield(searchOptions, {'InitialMeshSize', 'MaxMeshSize'});
defaultMeshOptions.MaxIterations = 1;
defaultMeshResult = pattern_search_cond4(eye(2), Q45, banks, ...
    @norm4_columns, normOptions, defaultMeshOptions);
assert(defaultMeshResult.options.InitialMeshSize == 0.05);
assert(defaultMeshResult.options.MaxMeshSize == 0.25);

% A Givens step may leave the canonical domain; return that accepted step.
oneStepOptions = struct('InitialMeshSize', pi/3, 'MaxMeshSize', pi/3, ...
    'MaxIterations', 1, 'Display', 'off', 'TrackHistory', true);
oneStepResult = pattern_search_cond4(eye(2), Q45, banks, ...
    @norm4_columns, normOptions, oneStepOptions);
canonicalFinal = canonicalize_rotations_by_trace_batch( ...
    oneStepResult.rawFinalRotations);
assert(oneStepResult.history.accepted(1,2));
assert(norm(oneStepResult.rawFinalRotations - canonicalFinal, 'fro') > 1e-3);
assert(isequal(oneStepResult.bestRotation, oneStepResult.rawFinalRotations));
assert(oneStepResult.bestValue == oneStepResult.rawFinalValues);
assert(strcmp(oneStepResult.terminationReasons{1}, 'MaxIterations'));
assert(oneStepResult.iterations(1) == 1);
assert(oneStepResult.functionCounts(1) > 1);

% The mesh threshold remains an independent stopping condition.
meshStopOptions = struct('InitialMeshSize', 1e-3, ...
    'MeshTolerance', 1e-2, 'MaxIterations', 3, 'Display', 'off');
meshStopResult = pattern_search_cond4(eye(2), Q45, banks, ...
    @norm4_columns, normOptions, meshStopOptions);
assert(strcmp(meshStopResult.terminationReasons{1}, 'MeshTolerance'));
assert(meshStopResult.iterations(1) == 1);

% A small but genuine improvement is accepted without a descent threshold.
smallAngle = 1e-9;
Qsmall = [cos(pi/8), -sin(pi/8); sin(pi/8), cos(pi/8)];
Qlower = [cos(pi/8 - smallAngle), -sin(pi/8 - smallAngle); ...
    sin(pi/8 - smallAngle), cos(pi/8 - smallAngle)];
smallInitial = compute_cond4_pq_fixed_starts(eye(2), Qsmall, banks, ...
    @norm4_columns, normOptions);
smallTrial = compute_cond4_pq_fixed_starts(eye(2), Qlower, banks, ...
    @norm4_columns, normOptions);
assert(smallTrial < smallInitial && smallInitial - smallTrial < 1e-8);
smallOptions = searchOptions;
smallOptions.InitialMeshSize = smallAngle;
smallOptions.MaxMeshSize = smallAngle;
smallOptions.MeshTolerance = 1e-12;
smallOptions.MaxIterations = 1;
smallResult = pattern_search_cond4(eye(2), Qsmall, banks, ...
    @norm4_columns, normOptions, smallOptions);
assert(smallResult.history.accepted(1,2));
assert(smallResult.rawFinalValues(1) < smallResult.initialValues(1));

% A progressive schedule must preserve the state of the surviving trajectory
% while pruning the other starts after the first stage.
Qprogressive = cat(3, Q45, eye(2), [cos(pi/8), -sin(pi/8); sin(pi/8), cos(pi/8)]);
progressiveOptions = struct('InitialMeshSize', pi/16, ...
    'MeshTolerance', 1e-8, 'MaxMeshSize', pi/4, ...
    'MeshExpansionFactor', 2, 'MeshContractionFactor', 0.5, ...
    'MaxIterations', 3, ...
    'Display', 'off', ...
    'StageIterationLimits', [1, 3], 'StageRetainedCounts', [3, 1]);
progressiveResult = pattern_search_cond4(P, Qprogressive, banks, ...
    @norm4_columns, normOptions, progressiveOptions);
assert(isequal(progressiveResult.finalRotations, ...
    progressiveResult.rawFinalRotations));
assert(isequal(progressiveResult.finalValues, progressiveResult.rawFinalValues));
[expectedBest, expectedIndex] = min(progressiveResult.rawFinalValues);
assert(progressiveResult.bestStartIndex == expectedIndex);
assert(progressiveResult.bestValue == expectedBest);
assert(isequal(progressiveResult.bestRotation, ...
    progressiveResult.rawFinalRotations(:,:,expectedIndex)));
directBest = compute_cond4_pq_fixed_starts(P, ...
    progressiveResult.bestRotation, banks, @norm4_columns, normOptions);
assert(abs(directBest - progressiveResult.bestValue) < 1e-12);
assert(sum(strcmp(progressiveResult.terminationReasons, 'StagePruned')) == 2);
assert(sum(progressiveResult.iterations == 3) == 1);
assert(all(progressiveResult.iterations( ...
    strcmp(progressiveResult.terminationReasons, 'StagePruned')) == 1));
assert(isequal(progressiveResult.options.StageIterationLimits, [1, 3]));
assert(isequal(progressiveResult.options.StageRetainedCounts, [3, 1]));
fullOptions = rmfield(progressiveOptions, ...
    {'StageIterationLimits', 'StageRetainedCounts'});
fullResult = pattern_search_cond4(P, Qprogressive, banks, ...
    @norm4_columns, normOptions, fullOptions);
assert(sum(progressiveResult.functionCounts) < sum(fullResult.functionCounts));
survivorIndex = find(~strcmp( ...
    progressiveResult.terminationReasons, 'StagePruned'));
assert(isscalar(survivorIndex));
assert(norm(progressiveResult.rawFinalRotations(:,:,survivorIndex) - ...
    fullResult.rawFinalRotations(:,:,survivorIndex), 'fro') < 1e-14);
assert(progressiveResult.rawFinalValues(survivorIndex) == ...
    fullResult.rawFinalValues(survivorIndex));
assert(progressiveResult.finalMeshSizes(survivorIndex) == ...
    fullResult.finalMeshSizes(survivorIndex));
assert(progressiveResult.functionCounts(survivorIndex) == ...
    fullResult.functionCounts(survivorIndex));

invalidStageOptions = progressiveOptions;
invalidStageOptions.StageRetainedCounts = [2, 1];
caughtInvalidSchedule = false;
try
    pattern_search_cond4(P, Qprogressive, banks, @norm4_columns, ...
        normOptions, invalidStageOptions);
catch exception
    caughtInvalidSchedule = strcmp(exception.identifier, ...
        'optimizer:pattern_search_cond4:InvalidStageSchedule');
end
assert(caughtInvalidSchedule);

disp('Pattern Search condition-number checks: OK');
