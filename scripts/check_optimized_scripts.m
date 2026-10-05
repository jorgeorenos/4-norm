% Check that example sweeps reuse fixed starts in batched power iterations.
oldVisibility = get(groot, 'DefaultFigureVisible');
set(groot, 'DefaultFigureVisible', 'off');
try
    profile clear;
    profile on;
    run(fullfile(fileparts(mfilename('fullpath')), 'optimizer_convergence.m'));
    profile off;
    info = profile('info');
    names = {info.FunctionTable.FunctionName};
    index = find(strcmp(names, 'power_norm4_multiple_starts'), 1);
    assert(~isempty(index) && info.FunctionTable(index).NumCalls == N_MATRICES);
    assert(isequal(size(norm_q_values), [N_MATRICES, numel(START_COUNTS)]));
    rng(SEED, 'twister');
    referenceQ = haar_so(DIMENSION);
    referenceState = rng;
    for j = 1:min(4, numel(START_COUNTS))
        referenceOptions = baseOptions;
        referenceOptions.NumRandomStarts = START_COUNTS(j);
        rng(referenceState);
        referenceValue = compute_induced_norm(referenceQ, @norm_4, referenceOptions);
        assert(abs(norm_q_values(1,j) - referenceValue) < 1e-12);
    end
    close all;

    profile clear;
    profile on;
    run(fullfile(fileparts(mfilename('fullpath')), 'optimerz_map.m'));
    profile off;
    info = profile('info');
    names = {info.FunctionTable.FunctionName};
    index = find(strcmp(names, 'power_norm4_multiple_starts'), 1);
    assert(~isempty(index) && ...
        info.FunctionTable(index).NumCalls == numel(iterationBudgets));
    assert(~any(strcmp(names, 'power_norm4_single_start')));
    assert(isequal(size(values), [size(X0,2), numel(iterationBudgets)]));
    for k = 1:numel(iterationBudgets)
        referenceOptions = options;
        referenceOptions.MaxIterations = iterationBudgets(k);
        for j = 1:size(X0,2)
            referenceValue = power_norm4_single_start(Q, @norm_4, ...
                X0(:,j), referenceOptions);
            assert(abs(values(j,k) - referenceValue) < 1e-12);
        end
    end
    close all;

    profile clear;
    profile on;
    run(fullfile(fileparts(mfilename('fullpath')), 'example_norm4_2x2.m'));
    profile off;
    info = profile('info');
    names = {info.FunctionTable.FunctionName};
    index = find(strcmp(names, 'norm_4'), 1);
    assert(~isempty(index) && info.FunctionTable(index).NumCalls < 500);
    assert(abs(qNorm4-angularReference) < 1e-8);
    close all;
catch err
    profile off;
    set(groot, 'DefaultFigureVisible', oldVisibility);
    rethrow(err);
end
set(groot, 'DefaultFigureVisible', oldVisibility);
disp('Optimized script checks: OK');
