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
    'MeshContractionFactor', 0.5, 'FunctionTolerance', 1e-10, ...
    'MaxIterations', 40, 'MaxFunctionEvaluations', 200, 'Display', 'off');
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

disp('Pattern Search condition-number checks: OK');
