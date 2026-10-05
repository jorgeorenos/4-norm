% Performance gate for the fixed SO(9), 1000-page, 100-start workload.
scriptDirectory = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDirectory);
addpath(fullfile(projectRoot, 'src', '4-norm'));
addpath(fullfile(projectRoot, 'src', 'helpers'));
rng(42, 'twister');
Q = zeros(9, 9, 1000);
for k = 1:1000
    Q(:,:,k) = haar_so(9);
end
state = rng;
options = struct('NumRandomStarts', 100, 'MaxIterations', 1000, ...
    'MaxWorkingMemoryMB', 128);
rng(state);
compute_induced_norm(Q, @norm_4, options); % Warm up
elapsed = zeros(1, 3);
for repeat = 1:3
    rng(state);
    timer = tic;
    values = compute_induced_norm(Q, @norm_4, options);
    elapsed(repeat) = toc(timer);
end
assert(isequal(size(values), [1000 1]) && all(isfinite(values)) && ...
    all(values >= 1-1e-12 & values <= 9^(1/4)+1e-12));
% A portable numerical check: the vector-only equivalent handle takes the
% same starts and should give the same estimates on a representative subset.
rng(state);
reference = compute_induced_norm(Q(:,:,1:20), @(x) norm_4(x), options);
rng(state);
fast = compute_induced_norm(Q(:,:,1:20), @norm_4, options);
assert(max(abs(reference - fast)) < 1e-12);
fprintf('100-start SO(9) performance median: %.6f s (%s)\n', ...
    median(elapsed), mat2str(elapsed));
