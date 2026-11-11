% Map kappa_4(P*Q) over 100000 canonical SO(3) Haar rotations.
% P is a reproducible Cholesky impact factor. The three panels reproduce
% one row of pairwise Euler-angle maps; each color is the minimum estimate
% within a projected angular bin over the omitted coordinate.
scriptDirectory = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptDirectory));
addpath(fullfile(projectRoot, 'src', '4-norm_handle'), '-begin');
addpath(fullfile(projectRoot, 'src', 'helpers'), '-begin');
addpath(fullfile(projectRoot, 'src', 'optimizer'), '-begin');

K = 3;
numSamples = 100000;
batchSize = 1000;
numRandomStarts = 75;
projectionLimitDegrees = 50;
projectionBinWidthDegrees = 0.5;
seeds = struct('Impact', 41, 'Haar', 42, 'Starts', 43);

% The fixed banks make every kappa_4(P*Q) estimate deterministic and keep
% its sampled directions identical for every canonical rotation.
normOptions = struct('MaxIterations', 1000, ...
    'NormTolerance', 1e-10, 'StationarityTolerance', 1e-8, ...
    'FeasibilityTolerance', 1e-12, 'MaxWorkingMemoryMB', 128, ...
    'ScreenIterations', 2, 'NumFinalists', 4);
normHandle = @norm4_columns;

rng(seeds.Impact, 'twister');
gaussianMatrix = randn(K);
SigmaE = gaussianMatrix*gaussianMatrix' + 0.25*eye(K);
P = chol((SigmaE + SigmaE')/2, 'lower');
banks = create_norm4_start_banks(K, numRandomStarts, seeds.Starts);

% The Haar seed is reset after constructing P. Consequently P, the start
% banks, the 100000 Q matrices, and the colors are reproducible on rerun.
rng(seeds.Haar, 'twister');
angles = zeros(3, numSamples);
conditionNumbers = zeros(numSamples, 1);
maximumOrthogonalityError = 0;
maximumDeterminantError = 0;

timer = tic;
for batchStart = 1:batchSize:numSamples
    batchEnd = min(batchStart + batchSize - 1, numSamples);
    batchIndices = batchStart:batchEnd;
    count = numel(batchIndices);
    rotations = zeros(K, K, count);
    for page = 1:count
        rotations(:,:,page) = haar_so(K);
    end

    % Right signed permutations leave the exact induced 4-norm condition
    % number invariant. Canonicalize before evaluating so one representative
    % is used consistently in the angle maps and in the fixed-bank estimate.
    canonicalRotations = canonicalize_so3_signed_permutations(rotations);
    conditionNumbers(batchIndices) = compute_cond4_pq_fixed_starts( ...
        P, canonicalRotations, banks, normHandle, normOptions);
    angles(:,batchIndices) = so3_zyx_euler_angles(canonicalRotations);

    for page = 1:count
        Q = canonicalRotations(:,:,page);
        maximumOrthogonalityError = max(maximumOrthogonalityError, ...
            norm(Q'*Q - eye(K), 'fro'));
        maximumDeterminantError = max(maximumDeterminantError, ...
            abs(det(Q) - 1));
    end
end
elapsedSeconds = toc(timer);

covarianceError = norm(P*P' - SigmaE, 'fro') / norm(SigmaE, 'fro');
assert(all(isfinite(conditionNumbers)) && all(conditionNumbers > 0));
assert(all(isfinite(angles(:))));
assert(maximumOrthogonalityError < 1e-12*max(1, K));
assert(maximumDeterminantError < 1e-12*max(1, K));
assert(covarianceError < 1e-12*max(1, K));

anglesDegrees = rad2deg(angles);
[minimumCondition, minimumIndex] = min(conditionNumbers);
alpha = anglesDegrees(1,:).';
beta = anglesDegrees(2,:).';
theta = anglesDegrees(3,:).';

alphaBetaMap = projected_minimum_map(alpha, beta, conditionNumbers, ...
    projectionLimitDegrees, projectionBinWidthDegrees);
alphaThetaMap = projected_minimum_map(alpha, theta, conditionNumbers, ...
    projectionLimitDegrees, projectionBinWidthDegrees);
betaThetaMap = projected_minimum_map(beta, theta, conditionNumbers, ...
    projectionLimitDegrees, projectionBinWidthDegrees);
logConditionValues = log10([alphaBetaMap.values; alphaThetaMap.values; ...
    betaThetaMap.values]);
colorLimits = [min(logConditionValues), max(logConditionValues)];
if diff(colorLimits) == 0
    colorLimits = colorLimits + [-0.01, 0.01];
end

fprintf('\nSO(3) Cholesky-factor objective map\n');
fprintf('Evaluated %d canonical Q matrices in %.3f seconds.\n', ...
    numSamples, elapsedSeconds);
fprintf('Fixed-bank kappa_4(P*Q) estimate range: [%.12f, %.12f].\n', ...
    minimumCondition, max(conditionNumbers));
fprintf('Sampled minimum: %.12f at [alpha beta theta] = [%.6f %.6f %.6f] degrees.\n', ...
    minimumCondition, alpha(minimumIndex), beta(minimumIndex), theta(minimumIndex));
fprintf('Maximum Q orthogonality error: %.3e | determinant error: %.3e.\n', ...
    maximumOrthogonalityError, maximumDeterminantError);
fprintf('Relative Cholesky covariance reconstruction error: %.3e.\n', ...
    covarianceError);

figureHandle = figure('Name', 'SO(3) Cholesky-factor condition map', 'Color', 'w');
layout = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
colormap(figureHandle, turbo(256));

axesHandles = gobjects(3, 1);
axesHandles(1) = nexttile;
draw_projected_map(axesHandles(1), alphaBetaMap, alpha(minimumIndex), ...
    beta(minimumIndex), '\alpha (degrees)', '\beta (degrees)', ...
    '\beta vs \alpha; min. over \theta', projectionLimitDegrees);
axesHandles(2) = nexttile;
draw_projected_map(axesHandles(2), alphaThetaMap, alpha(minimumIndex), ...
    theta(minimumIndex), '\alpha (degrees)', '\theta (degrees)', ...
    '\theta vs \alpha; min. over \beta', projectionLimitDegrees);
axesHandles(3) = nexttile;
draw_projected_map(axesHandles(3), betaThetaMap, beta(minimumIndex), ...
    theta(minimumIndex), '\beta (degrees)', '\theta (degrees)', ...
    '\theta vs \beta; min. over \alpha', projectionLimitDegrees);
for axesIndex = 1:numel(axesHandles)
    clim(axesHandles(axesIndex), colorLimits);
end
colorbarHandle = colorbar(axesHandles(3));
colorbarHandle.Layout.Tile = 'east';
ylabel(colorbarHandle, 'log_{10}(minimum estimated \kappa_4(PQ))');
sgtitle(layout, sprintf(['Cholesky P | %d canonical Q matrices | ', ...
    'shared color scale'], numSamples));

function map = projected_minimum_map(x, y, values, limitDegrees, binWidthDegrees)
    edges = -limitDegrees:binWidthDegrees:limitDegrees;
    centers = edges(1:end-1) + binWidthDegrees/2;
    binCount = numel(centers);
    xBin = discretize(x, edges);
    yBin = discretize(y, edges);
    valid = ~isnan(xBin) & ~isnan(yBin);
    linearBins = sub2ind([binCount, binCount], yBin(valid), xBin(valid));
    minimumValues = accumarray(linearBins, values(valid), ...
        [binCount*binCount, 1], @min, NaN);
    [row, column] = find(reshape(isfinite(minimumValues), binCount, binCount));
    map = struct('x', centers(column).', 'y', centers(row).', ...
        'values', minimumValues(sub2ind([binCount, binCount], row, column)));
end

function draw_projected_map(axesHandle, map, minimumX, minimumY, xLabel, yLabel, titleText, limitDegrees)
    scatter(axesHandle, map.x, map.y, 7, log10(map.values), 'filled', ...
        'MarkerEdgeColor', 'none');
    hold(axesHandle, 'on');
    plot(axesHandle, minimumX, minimumY, 'r*', 'MarkerSize', 9, 'LineWidth', 1.2);
    hold(axesHandle, 'off');
    axis(axesHandle, 'square');
    xlim(axesHandle, [-limitDegrees, limitDegrees]);
    ylim(axesHandle, [-limitDegrees, limitDegrees]);
    grid(axesHandle, 'on');
    xlabel(axesHandle, xLabel);
    ylabel(axesHandle, yLabel);
    title(axesHandle, titleText);
end