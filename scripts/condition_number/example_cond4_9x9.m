% Compare the induced 4-norm condition number estimate for 1000 random
% matrices in SO(9). Each field of cond4_numbers contains one estimate per
% matrix for a fixed number of random starts. The same Q is used for every
% start count. For 50 or more starts, screen each for up to two iterations
% and refine the best four. Smaller groups retain full iteration to avoid
% screening bias.

scriptDirectory = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptDirectory));
addpath(fullfile(projectRoot, 'src', '4-norm'));
addpath(fullfile(projectRoot, 'src', 'helpers'));

rng(42, 'twister');

n = 9;
numMatrices = 1000000;
startCounts = [10, 25, 50, 75,100, 250, 500, 750, 1000];
baseOptions = struct('MaxIterations', 1000, 'MaxWorkingMemoryMB', 128);
screenIterations = 3;
numFinalists = 3;
screenMinStarts = 50;

% Each field is a numMatrices-by-1 vector of condition number estimates.
cond4_numbers = struct();
fieldNames = cell(1, numel(startCounts));
elapsedByStartCount = zeros(1, numel(startCounts));

fprintf('Computing condition numbers for %d matrices in SO(%d)...\n', ...
    numMatrices, n);
totalTimer = tic;

rotations = zeros(n, n, numMatrices);
for k = 1:numMatrices
    rotations(:,:,k) = haar_so(n);
end
for j = 1:numel(startCounts)
    fieldNames{j} = sprintf('starts%d', startCounts(j));
    options = baseOptions;
    options.NumRandomStarts = startCounts(j);
    if startCounts(j) >= screenMinStarts
        options.ScreenIterations = screenIterations;
        options.NumFinalists = numFinalists;
    end
    startTimer = tic;
    cond4_numbers.(fieldNames{j}) = compute_4_cond(rotations, options);
    elapsedByStartCount(j) = toc(startTimer);
    fprintf('  Completed batch with %d random starts\n', startCounts(j));
end

totalElapsedTime = toc(totalTimer);
theoreticalUpperBound = sqrt(n);

fprintf('\n--- Results by number of random starts ---\n');
referenceValues = cond4_numbers.(fieldNames{1});
means = zeros(size(startCounts));
medians = zeros(size(startCounts));
standardDeviations = zeros(size(startCounts));
for j = 1:numel(startCounts)
    values = cond4_numbers.(fieldNames{j});
    means(j) = mean(values);
    medians(j) = median(values);
    standardDeviations(j) = std(values);
    maxDifference = max(abs(values - referenceValues));
    fprintf(['%4d random starts: mean %.12f | median %.12f | ' ...
        'std %.12f | %.3f s total | max. difference from starts1 %.3e\n'], ...
        startCounts(j), means(j), medians(j), ...
        standardDeviations(j), elapsedByStartCount(j), maxDifference);
    assert(all(values >= 1 - 1e-9 & ...
        values <= theoreticalUpperBound + 1e-9), ...
        'An estimate falls outside the theoretical bounds.');
end
fprintf('Total elapsed time: %.3f seconds\n', totalElapsedTime);

% Common limits and bins make the nine empirical distributions comparable.
binEdges = linspace(1 - 1e-9, theoreticalUpperBound + 1e-9, 31);
figure('Name', 'Condition number estimates in SO(9)');
tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for j = 1:numel(startCounts)
    values = cond4_numbers.(fieldNames{j});
    nexttile;
    histogram(values, binEdges, 'FaceColor', [0.85, 0.40, 0.10], ...
        'EdgeColor', [0.25, 0.25, 0.25]);
    grid on;
    xlim([binEdges(1), binEdges(end)]);
    xlabel('Estimated condition number');
    ylabel('Number of matrices');
    title(sprintf(['%d random starts\nmean = %.6f | median = %.6f\n' ...
        'standard deviation = %.6f | time = %.1f s'], startCounts(j), ...
        means(j), medians(j), standardDeviations(j), ...
        elapsedByStartCount(j)));
end

sgtitle(sprintf(['Condition number estimates for %d matrices in SO(%d)\n' ...
    'Same matrices; screened estimates for %d or more starts'], ...
    numMatrices, n, screenMinStarts));
