% Reproducible timings of the induced 4-norm; matrix generation is excluded.
scriptDirectory = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDirectory);
addpath(fullfile(projectRoot, 'src', '4-norm'));
addpath(fullfile(projectRoot, 'src', 'helpers'));
rng(42, 'twister');

n = 9;
numMatrices = 1000;
numRepeats = 3;
options = struct('NumRandomStarts', 100, 'MaxIterations', 1000, ...
    'MaxWorkingMemoryMB', 128);
Qbatch = zeros(n,n,numMatrices);
for k = 1:numMatrices
    Qbatch(:,:,k) = haar_so(n);
end
initialState = rng;
compute_induced_norm(Qbatch, @norm_4, options); % Warm up MATLAB.
elapsed = zeros(numRepeats,1);
for repeat = 1:numRepeats
    rng(initialState);
    timer = tic;
    values = compute_induced_norm(Qbatch, @norm_4, options);
    elapsed(repeat) = toc(timer);
end
assert(isequal(size(values), [numMatrices 1]) && all(isfinite(values)));
fprintf('MATLAB %s\n', version);
fprintf('%d matrices in SO(%d), %d starts: median %.6f s\n', ...
    numMatrices, n, options.NumRandomStarts, median(elapsed));
disp(elapsed);
