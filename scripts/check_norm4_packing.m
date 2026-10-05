% Packed and scalar trajectories must retain each original start's result.
scriptDirectory = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDirectory);
addpath(fullfile(projectRoot, 'src', '4-norm'));
addpath(fullfile(projectRoot, 'src', 'helpers'));
rng(71, 'twister');
n = 9;
numMatrices = 8;
numStarts = 100;
Q = zeros(n, n, numMatrices);
for page = 1:numMatrices
    Q(:,:,page) = haar_so(n);
end
starts = randn(n, numStarts, numMatrices);
options = struct('MaxIterations', 150, 'NormTolerance', 1e-10, ...
    'StationarityTolerance', 1e-8, 'FeasibilityTolerance', 1e-12);
packed = power_norm4_multiple_starts(Q, @norm_4, starts, options);
assert(isequal(size(packed), [numStarts numMatrices]));
for page = 1:numMatrices
    for start = 1:numStarts
        scalar = power_norm4_multiple_starts(Q(:,:,page), @norm_4, ...
            starts(:,start,page), options);
        assert(abs(packed(start,page) - scalar) < 1e-12);
    end
end
% The general vector-only handle path must remain supported.
vectorOnly = power_norm4_multiple_starts(Q, @(x) norm_4(x), starts, options);
assert(max(abs(packed(:)-vectorOnly(:))) < 1e-12);
disp('Packed 4-norm trajectories: OK');
