% Basic checks of the induced 4-norm estimate (base MATLAB).
scriptDirectory = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDirectory);
addpath(fullfile(projectRoot, 'src', '4-norm'));
addpath(fullfile(projectRoot, 'src', 'helpers'));

assert(abs(norm_4([1; 2]) - 17^(1/4)) < 1e-14);
assert(norm_4([0, 0]) == 0);
assert(norm_4([1, 2]) == norm_4([1; 2]));
assert(isfinite(norm_4([1e200; 1e200])));
assert(norm_4([1e-200; 1e-200]) > 0);

opts = struct('NumRandomStarts', 3, 'MaxIterations', 100);
for n = [1 2 5]
    [value, x, info] = compute_induced_norm(eye(n), @norm_4, opts);
    assert(abs(value - 1) < 1e-12);
    assert(abs(norm_4(x) - 1) < 1e-12 && info.bestConverged);
end
Q = [1, -1; 1, 1] / sqrt(2);
[value, x, info] = compute_induced_norm(Q, @norm_4, opts);
assert(abs(value - 2^(1/4)) < 1e-8);
assert(abs(value - norm_4(Q*x)) < 1e-12 && info.bestConverged);
P = [0, -1; 1, 0];
assert(abs(compute_induced_norm(P, @norm_4, opts) - 1) < 1e-12);
for n = [2 5 10]
    Qhaar = haar_so(n);
    [estimate, x, diagnostics] = compute_induced_norm(Qhaar, @norm_4, opts);
    assert(isfinite(estimate) && abs(estimate - norm_4(Qhaar*x)) < 1e-12);
    assert(abs(norm_4(x) - 1) < 1e-12);
    assert(estimate >= 1-1e-12 && estimate <= n^(1/4)+1e-12);
    assert(diagnostics.numStarts == n + opts.NumRandomStarts);
    assert(isfinite(diagnostics.bestStationarityResidual));
end
rng(42,'twister');
Q1 = haar_so(2);
[a, xa] = compute_induced_norm(Q1, @norm_4, opts);
rng(42,'twister');
Q2 = haar_so(2);
[b, xb] = compute_induced_norm(Q2, @norm_4, opts);
assert(isequal(Q1,Q2) && isequal(a,b) && isequal(xa,xb));
rng(42,'twister');
Q3 = haar_so(2);
[c, xc] = compute_induced_norm(Q3, @(x) norm_4(x), opts);
assert(isequal(Q1,Q3) && abs(a-c) < 1e-12 && norm_4(xa-xc) < 1e-12);
historyOptions = struct('NumRandomStarts', 2, 'MaxIterations', 100, ...
    'StoreHistory', true);
[~, ~, historyInfo] = compute_induced_norm(eye(2), @norm_4, historyOptions);
assert(numel(historyInfo.histories) == historyInfo.numStarts && ...
    all(cellfun(@(h) ~isempty(h), historyInfo.histories)));
stateBefore = rng;
compute_induced_norm(Q1, @norm_4, opts);
stateAfter = rng;
assert(~isequal(stateBefore.State, stateAfter.State));
short = struct('NumRandomStarts', 0, 'MaxIterations', 1, ...
    'NormTolerance', 1e-30, 'StationarityTolerance', 1e-30);
lastwarn('');
[~, ~, limited] = compute_induced_norm(Q, @norm_4, short);
[~, warningId] = lastwarn;
assert(~limited.bestConverged && strcmp(limited.bestTerminationReason,'maxIterations'));
assert(strcmp(warningId,'compute_induced_norm:NotConverged'));
assert_error(@() norm_4(ones(2,2)), 'norm_4:InvalidVector');
assert_error(@() norm_4([1;Inf]), 'norm_4:InvalidVector');
assert_error(@() compute_induced_norm(ones(2), @norm_4), 'compute_induced_norm:NotSO');
assert_error(@() compute_induced_norm(eye(2), @norm_4, struct('Bad',1)), ...
    'compute_induced_norm:UnknownOption');
assert_error(@() compute_induced_norm(eye(2), @norm_4, ...
    struct('NumRandomStarts',-1)), 'compute_induced_norm:InvalidCount');
assert_error(@() compute_induced_norm(eye(2), @norm_4, ...
    struct('NormTolerance',0)), 'compute_induced_norm:InvalidTolerance');
assert_error(@() compute_induced_norm(eye(2), @(x) 0, opts), ...
    'compute_induced_norm:InvalidHandleOutput');
assert_error(@() compute_induced_norm(eye(2), @(x) NaN, opts), ...
    'compute_induced_norm:InvalidHandleOutput');
disp('Basic 4-norm checks: OK');

function assert_error(thunk, identifier)
    try
        thunk();
    catch err
        assert(strcmp(err.identifier, identifier), err.message);
        return;
    end
    error('check_norm4:ExpectedError', 'Expected error %s.', identifier);
end
