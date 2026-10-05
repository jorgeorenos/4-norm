% Time the screened induced 4-norm estimate for 1000 matrices and 100 starts.
% Every start receives up to four iterations; the best three per matrix
% are then refined. Screening can miss a late-improving start.
scriptDirectory = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDirectory);
addpath(fullfile(projectRoot, 'src', '4-norm'));
addpath(fullfile(projectRoot, 'src', 'helpers'));

rng(42, 'twister');
n = 9;
numMatrices = 1000;
numStarts = 100;
rotations = zeros(n, n, numMatrices);
for k = 1:numMatrices
    rotations(:,:,k) = haar_so(n);
end

options = struct('NumRandomStarts', numStarts, 'MaxIterations', 1000, ...
    'MaxWorkingMemoryMB', 128, 'ScreenIterations', 2, 'NumFinalists', 4);
timer = tic;
estimates = compute_4_norm(rotations, options);
elapsedSeconds = toc(timer);

fprintf('Matrices Q: %d\n', numMatrices);
fprintf('Random starting vectors per matrix: %d\n', numStarts);
fprintf('Computation time: %.6f seconds\n', elapsedSeconds);
