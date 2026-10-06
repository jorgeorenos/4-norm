% Time the screened condition number estimate for 1000 matrices and 100
% starts. Every start receives up to two iterations; the best four per
% matrix are then refined. Screening can miss a late-improving start.
scriptDirectory = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptDirectory));
addpath(fullfile(projectRoot, 'src', '4-norm'));
addpath(fullfile(projectRoot, 'src', 'helpers'));

rng(42, 'twister');
n = 9;
numMatrices = 100000;
numStarts = 75;
rotations = zeros(n, n, numMatrices);
for k = 1:numMatrices
    rotations(:,:,k) = haar_so(n);
end

options = struct('NumRandomStarts', numStarts, 'MaxIterations', 1000, ...
    'MaxWorkingMemoryMB', 128, 'ScreenIterations', 3, 'NumFinalists', 3);
timer = tic;
estimates = compute_4_cond(rotations, options);
elapsedSeconds = toc(timer);

fprintf('Matrices Q: %d\n', numMatrices);
fprintf('Random starting vectors per matrix: %d\n', numStarts);
fprintf('Computation time: %.6f seconds\n', elapsedSeconds);
