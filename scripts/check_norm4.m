% Basic checks of the induced 4-norm estimate (base MATLAB).
scriptDirectory = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDirectory);
addpath(fullfile(projectRoot, 'src', '4-norm'));
addpath(fullfile(projectRoot, 'src', 'helpers'));

assert(abs(norm4_columns([1; 2]) - 17^(1/4)) < 1e-14);
assert(all(norm4_columns([0, 0]) == 0));
assert(isfinite(norm4_columns([1e200; 1e200])));
assert(norm4_columns([1e-200; 1e-200]) > 0);
normSamples = cat(3, [0, 1e200; 0, -1e200], [1e-200, 3; -1e-200, 4]);
normValues = norm4_columns(normSamples);
assert(isequal(size(normValues), [1 2 2]));
for page = 1:2
    for column = 1:2
        assert(normValues(1,column,page) == norm4_columns(normSamples(:,column,page)));
    end
end
assert(isequal(norm4_columns([1 2]), [1 2]));

rng(42, 'twister');
opts = struct('NumRandomStarts', 20, 'MaxIterations', 100);
for n = [1 2 5 10]
    assert(abs(compute_4_norm(eye(n), opts) - 1) < 1e-12);
    directions = randn(n, 25);
    X = directions ./ norm4_columns(directions);
    assert(isequal(size(X), [n 25]) && all(isfinite(X(:))));
    assert(max(abs(sum(abs(X).^4,1).^(1/4)-1)) < 1e-12);
    Qhaar = haar_so(n);
    assert(norm(Qhaar'*Qhaar-eye(n),'fro') < 1e-12*max(1,n));
    assert(abs(det(Qhaar)-1) < 1e-12*max(1,n));
    estimate = compute_4_norm(Qhaar, opts);
    assert(isfinite(estimate) && estimate >= 1-1e-12 && estimate <= n^(1/4)+1e-12);
end
Q = [1, -1; 1, 1] / sqrt(2);
P = [0, -1; 1, 0];
batch = cat(3, eye(2), Q, P, haar_so(2));
startsState = rng;
values = compute_4_norm(batch, opts);
batchEndState = rng;
assert(isequal(size(values), [4 1]));
assert(max(abs(values(1:3)-[1; 2^(1/4); 1])) < 1e-8);
rng(startsState);
sequential = zeros(4,1);
for k = 1:4
    sequential(k) = compute_4_norm(batch(:,:,k), opts);
end
assert(isequal(values, sequential) && isequal(rng, batchEndState));
rng(startsState);
assert(isequal(values, compute_4_norm(batch, opts)));
% Page removal must preserve ordering and random streams across block sizes.
for n = [1 2 5 10]
    matrices = zeros(n,n,7);
    matrices(:,:,1) = eye(n);
    for page = 2:7
        matrices(:,:,page) = haar_so(n);
    end
    initialState = rng;
    batched = compute_4_norm(matrices, opts);
    finalState = rng;
    tinyBlocks = opts;
    tinyBlocks.MaxWorkingMemoryMB = 0.001; % Force one matrix per block.
    rng(initialState);
    blocked = compute_4_norm(matrices, tinyBlocks);
    assert(max(abs(batched-blocked)) < 1e-12 && isequal(rng,finalState));
    rng(initialState);
    for page = 1:7
        rawStarts = randn(n,opts.NumRandomStarts);
        singleValues = zeros(1,opts.NumRandomStarts);
        for column = 1:opts.NumRandomStarts
            singleValues(column) = power_norm4_multiple_starts( ...
                matrices(:,:,page), rawStarts(:,column), ...
                struct('MaxIterations', opts.MaxIterations, 'NormTolerance', 1e-10, ...
                'StationarityTolerance', 1e-8, 'FeasibilityTolerance', 1e-12));
        end
        assert(abs(batched(page)-max(singleValues)) < 1e-12);
    end
    assert(isequal(rng,finalState));
end
assert(~isequal(startsState.State, batchEndState.State));
assert(isequal(compute_4_norm(ones(1,1,3), opts), ones(3,1)));
short = struct('NumRandomStarts', 1, 'MaxIterations', 1);
assert(isfinite(compute_4_norm(Q, short)));
assert(nargout('compute_4_norm') == 1);
assert(nargout('power_norm4_multiple_starts') == 1);
assert(isempty(which('compute_induced_norm')));
assert(isempty(which('norm_4')));
assert(isempty(which('generate_unit_l4_vectors')));
assert(isempty(which('power_norm4_single_start')));
assert_error(@() compute_4_norm(zeros(2,3,4)), 'compute_4_norm:InvalidMatrix');
assert_error(@() compute_4_norm(zeros(2,2,2,2)), 'compute_4_norm:InvalidMatrix');
assert_error(@() compute_4_norm(cat(3, eye(2), ones(2))), 'compute_4_norm:NotSO');
assert_error(@() compute_4_norm(eye(2), struct('StoreHistory',true)), ...
    'compute_4_norm:UnknownOption');
assert_error(@() compute_4_norm(eye(2), ...
    struct('MaxWorkingMemoryMB',0)), 'compute_4_norm:InvalidMemoryBudget');
assert_error(@() compute_4_norm(ones(2)), 'compute_4_norm:NotSO');
assert_error(@() compute_4_norm(eye(2), struct('Bad',1)), ...
    'compute_4_norm:UnknownOption');
assert_error(@() compute_4_norm(eye(2), ...
    struct('NumRandomStarts',-1)), 'compute_4_norm:InvalidCount');
assert_error(@() compute_4_norm(eye(2), ...
    struct('NumRandomStarts',0)), 'compute_4_norm:InvalidCount');
assert_error(@() compute_4_norm(eye(2), ...
    struct('NormTolerance',0)), 'compute_4_norm:InvalidTolerance');
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
