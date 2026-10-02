% Batch estimation of the induced 4-norm for 1000 random matrices in SO(9).
% This script measures the total wall-clock time required to run the
% implemented power-method estimator on a reproducible ensemble.

scriptDirectory = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDirectory);
addpath(fullfile(projectRoot, 'src', '4-norm'));
addpath(fullfile(projectRoot, 'src', 'helpers'));

rng(42, 'twister');

n = 9;
numMatrices = 1000;
normHandle = @norm_4;
options = struct('NumRandomStarts', 1000, 'MaxIterations', 1000);

normValues = zeros(numMatrices, 1);

fprintf('Computing induced 4-norm for %d matrices in SO(%d)...\n', numMatrices, n);
tic;

for k = 1:numMatrices
    Q = haar_so(n);
    qNorm4 = compute_induced_norm(Q, normHandle, options);
    normValues(k) = qNorm4;
    if mod(k, 100) == 0
        fprintf('  Completed %d / %d\n', k, numMatrices);
    end
end

elapsedTime = toc;

fprintf('\n--- Results ---\n');
fprintf('Total elapsed time: %.3f seconds (%.4f s per matrix)\n', ...
    elapsedTime, elapsedTime / numMatrices);
fprintf('Minimum ||Q||_4: %.12f\n', min(normValues));
fprintf('Maximum ||Q||_4: %.12f\n', max(normValues));
fprintf('Mean ||Q||_4:   %.12f\n', mean(normValues));
fprintf('Median ||Q||_4: %.12f\n', median(normValues));
fprintf('Theoretical bounds: [1, %.12f]\n', n^(1/4));
fprintf('All within bounds: %d\n', ...
    all(normValues >= 1 - 1e-12 & normValues <= n^(1/4) + 1e-12));
