% Basic checks of the induced 4-norm condition number estimate (base MATLAB).
scriptDirectory = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptDirectory));
addpath(fullfile(projectRoot, 'src', '4-norm'));
addpath(fullfile(projectRoot, 'src', 'helpers'));

rng(42, 'twister');
opts = struct('NumRandomStarts', 40, 'MaxIterations', 200);

% Exact values in low dimensions: kappa_4(I) = 1, kappa_4 of the 45-degree
% rotation is sqrt(2) (both factors equal 2^(1/4)), and the signed
% permutation has kappa_4 = 1.
assert(abs(compute_4_cond(1, opts) - 1) < 1e-12);
assert(abs(compute_4_cond(eye(2), opts) - 1) < 1e-10);
Q45 = [1, -1; 1, 1] / sqrt(2);
assert(abs(compute_4_cond(Q45, opts) - sqrt(2)) < 1e-6);
P = [0, -1; 1, 0];
assert(abs(compute_4_cond(P, opts) - 1) < 1e-10);

% Independent angular reference for a fixed rotation in dimension 2.
angle = 0.37;
Qref = [cos(angle), -sin(angle); sin(angle), cos(angle)];
theta = (0:99999)*(2*pi/100000);
directions = [cos(theta); sin(theta)];
X = directions ./ norm4_columns(directions);
referenceKappa = max(norm4_columns(Qref*X)) * max(norm4_columns(Qref'*X));
assert(abs(compute_4_cond(Qref, opts) - referenceKappa) < 1e-5);

% Haar batches: theoretical bounds 1 <= kappa_4(Q) <= n^(1/2) in SO(n),
% reproducibility, and equality with sequential per-matrix calls and with
% an explicit interleaved batch of Q and its transpose.
for n = [1 2 5 10]
    matrices = zeros(n, n, 7);
    for page = 1:7
        matrices(:,:,page) = haar_so(n);
    end
    initialState = rng;
    kappas = compute_4_cond(matrices, opts);
    finalState = rng;
    assert(isequal(size(kappas), [7 1]) && all(isfinite(kappas)));
    assert(all(kappas >= 1 - 1e-9 & kappas <= sqrt(n) + 1e-9));
    rng(initialState);
    assert(isequal(kappas, compute_4_cond(matrices, opts)));
    assert(isequal(rng, finalState));
    rng(initialState);
    sequential = zeros(7, 1);
    for page = 1:7
        sequential(page) = compute_4_cond(matrices(:,:,page), opts);
    end
    assert(isequal(kappas, sequential) && isequal(rng, finalState));
    rng(initialState);
    interleaved = zeros(n, n, 14);
    interleaved(:,:,1:2:end) = matrices;
    interleaved(:,:,2:2:end) = permute(matrices, [2 1 3]);
    values = compute_4_norm(interleaved, opts);
    assert(isequal(kappas, values(1:2:end) .* values(2:2:end)));
    assert(isequal(rng, finalState));
end

% A scalar input returns a scalar; one page returns a one-element vector.
assert(isscalar(compute_4_cond(Q45, opts)));
assert(isequal(size(compute_4_cond(cat(3, Q45, P), opts)), [2 1]));
assert(nargout('compute_4_cond') == 1);

assert_error(@() compute_4_cond(zeros(2,3,4)), 'compute_4_cond:InvalidMatrix');
assert_error(@() compute_4_cond(zeros(2,2,2,2)), 'compute_4_cond:InvalidMatrix');
assert_error(@() compute_4_cond([]), 'compute_4_cond:InvalidMatrix');
assert_error(@() compute_4_cond(single(eye(2))), 'compute_4_cond:InvalidMatrix');
assert_error(@() compute_4_cond(eye(2), struct('Bad',1)), ...
    'compute_4_norm:UnknownOption');
assert_error(@() compute_4_cond(eye(2), ...
    struct('NumRandomStarts',0)), 'compute_4_norm:InvalidCount');
disp('Basic condition number checks: OK');

function assert_error(thunk, identifier)
    try
        thunk();
    catch err
        assert(strcmp(err.identifier, identifier), err.message);
        return;
    end
    error('check_cond4:ExpectedError', 'Expected error %s.', identifier);
end
