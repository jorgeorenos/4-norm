% Compare the induced 4-norm estimate for 1000 random matrices in SO(9).
% Each field of norm4_norms contains one estimate per matrix for a fixed
% number of random starts. The same Q is used for every start count.

scriptDirectory = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDirectory);
addpath(fullfile(projectRoot, 'src', '4-norm'));
addpath(fullfile(projectRoot, 'src', 'helpers'));

rng(42, 'twister');

n = 9;
numMatrices = 1000;
startCounts = [1, 10, 25, 50, 100, 250, 500, 750, 1000];
normHandle = @norm_4;
baseOptions = struct('MaxIterations', 1000);

% MATLAB identifiers cannot start with a number, hence norm4_norms rather
% than 4_norms. Each field is a numMatrices-by-1 vector.
norm4_norms = struct();
fieldNames = cell(1, numel(startCounts));
elapsedByStartCount = zeros(1, numel(startCounts));
for j = 1:numel(startCounts)
    fieldNames{j} = sprintf('starts%d', startCounts(j));
    norm4_norms.(fieldNames{j}) = zeros(numMatrices, 1);
end

fprintf('Computing induced 4-norms for %d matrices in SO(%d)...\n', ...
    numMatrices, n);
totalTimer = tic;

for k = 1:numMatrices
    Q = haar_so(n);
    for j = 1:numel(startCounts)
        options = baseOptions;
        options.NumRandomStarts = startCounts(j);
        startTimer = tic;
        norm4_norms.(fieldNames{j})(k) = ...
            compute_induced_norm(Q, normHandle, options);
        elapsedByStartCount(j) = elapsedByStartCount(j) + toc(startTimer);
    end
    if mod(k, 100) == 0
        fprintf('  Completed %d / %d matrices\n', k, numMatrices);
    end
end

totalElapsedTime = toc(totalTimer);
theoreticalUpperBound = n^(1/4);

fprintf('\n--- Results by number of random starts ---\n');
referenceValues = norm4_norms.(fieldNames{1});
for j = 1:numel(startCounts)
    values = norm4_norms.(fieldNames{j});
    maxDifference = max(abs(values - referenceValues));
    fprintf(['%4d random starts: mean %.12f | median %.12f | ' ...
        'std %.12f | %.3f s total | max. difference from starts1 %.3e\n'], ...
        startCounts(j), mean(values), median(values), ...
        std(values), elapsedByStartCount(j), maxDifference);
    assert(all(values >= 1 - 1e-12 & ...
        values <= theoreticalUpperBound + 1e-12), ...
        'An estimate falls outside the theoretical bounds.');
end
fprintf('Total elapsed time: %.3f seconds\n', totalElapsedTime);

% Common limits and bins make the six empirical distributions comparable.
binEdges = linspace(1 - 1e-12, theoreticalUpperBound + 1e-12, 31);
figure('Name', 'Induced 4-norm estimates in SO(9)');
tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for j = 1:numel(startCounts)
    values = norm4_norms.(fieldNames{j});
    nexttile;
    histogram(values, binEdges, 'FaceColor', [0.85, 0.40, 0.10], ...
        'EdgeColor', [0.25, 0.25, 0.25]);
    grid on;
    xlim([binEdges(1), binEdges(end)]);
    xlabel('Estimated induced 4-norm');
    ylabel('Number of matrices');
    title(sprintf(['%d random starts\nmean = %.6f | median = %.6f\n' ...
        'standard deviation = %.6f | time = %.1f s'], startCounts(j), ...
        mean(values), median(values), std(values), ...
        elapsedByStartCount(j)));
end

sgtitle(sprintf(['Induced 4-norm estimates for %d matrices in SO(%d)\n' ...
    'Same matrices, varying number of random starts'], numMatrices, n));
