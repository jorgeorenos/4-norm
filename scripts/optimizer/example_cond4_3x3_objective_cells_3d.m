% Visualize binned kappa_4(P*Q) estimates in canonical SO(3) Euler space.
% Each visible dot is one occupied [alpha, beta, theta] angular cell and
% its color is the minimum fixed-bank estimate found in that cell.
scriptDirectory = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptDirectory));
addpath(fullfile(projectRoot, 'src', '4-norm_handle'), '-begin');
addpath(fullfile(projectRoot, 'src', 'helpers'), '-begin');
addpath(fullfile(projectRoot, 'src', 'optimizer'), '-begin');

K = 3;
numSamples = 10000;
batchSize = 1000;
numRandomStarts = 75;
angleLimitDegrees = 50;
cellWidthDegrees = 1.755;
maximumDisplayedCells = 15000;
seeds = struct('Impact', 41, 'Haar', 42, 'Starts', 43, 'Display', 44);

% Fixed start banks turn the multistart estimate into the same deterministic
% numerical objective for every Q. Screening remains an explicit heuristic.
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
alpha = anglesDegrees(1,:).';
beta = anglesDegrees(2,:).';
theta = anglesDegrees(3,:).';
[minimumCondition, minimumIndex] = min(conditionNumbers);

cellMap = occupied_condition_cells(alpha, beta, theta, conditionNumbers, ...
    angleLimitDegrees, cellWidthDegrees);
occupiedCount = numel(cellMap.values);
displayedCount = min(maximumDisplayedCells, occupiedCount);
displayStream = RandStream('mt19937ar', 'Seed', seeds.Display);
displayPositions = randperm(displayStream, occupiedCount, displayedCount);

figureHandle = figure('Name', 'SO(3) celdas ocupadas de condicion de Cholesky', 'Color', 'w');
axesHandle = axes(figureHandle);
scatter3(axesHandle, cellMap.alpha(displayPositions), ...
    cellMap.beta(displayPositions), cellMap.theta(displayPositions), 13, ...
    log10(cellMap.values(displayPositions)), 'filled', 'MarkerEdgeColor', 'none');
hold(axesHandle, 'on');
plot3(axesHandle, alpha(minimumIndex), beta(minimumIndex), theta(minimumIndex), ...
    'rp', 'MarkerSize', 11, 'MarkerFaceColor', 'r', 'LineWidth', 1.0);
hold(axesHandle, 'off');
colormap(figureHandle, turbo(256));
colorLimits = log10([min(cellMap.values), max(cellMap.values)]);
if diff(colorLimits) == 0
    colorLimits = colorLimits + [-0.01, 0.01];
end
clim(axesHandle, colorLimits);
colorbarHandle = colorbar(axesHandle);
ylabel(colorbarHandle, 'log_{10}(minimo estimado de \kappa_4(PQ))');
annotation(figureHandle, 'textbox', [0.77, 0.91, 0.10, 0.06], ...
    'String', {'Celdas ocupadas', ...
    sprintf('Mejor \\kappa = %.6f', minimumCondition)}, 'Interpreter', 'tex', ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'middle', 'FontSize', 9);

axis(axesHandle, 'vis3d');
xlim(axesHandle, [-angleLimitDegrees, angleLimitDegrees]);
ylim(axesHandle, [-angleLimitDegrees, angleLimitDegrees]);
zlim(axesHandle, [-angleLimitDegrees, angleLimitDegrees]);
grid(axesHandle, 'on');
box(axesHandle, 'on');
view(axesHandle, [-38, 22]);
xlabel(axesHandle, '\alpha (grados)');
ylabel(axesHandle, '\beta (grados)');
zlabel(axesHandle, '\theta (grados)');
annotation(figureHandle, 'textbox', [0.02, 0.96, 0.78, 0.03], ...
    'String', sprintf(['K=%d · %s rotaciones · %s/%s celdas ocupadas mostradas'], ...
    K, thousands_label(numSamples), thousands_label(displayedCount), ...
    thousands_label(occupiedCount)), 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');

fprintf('\nMapa de celdas ocupadas SO(3) con factor de Cholesky\n');
fprintf('Se evaluaron %d matrices Q canonicas en %.3f segundos.\n', ...
    numSamples, elapsedSeconds);
fprintf('Se muestran %d de %d celdas tridimensionales ocupadas.\n', ...
    displayedCount, occupiedCount);
fprintf('Minimo muestral fijo de kappa_4(P*Q): %.12f en [alpha beta theta] = [%.6f %.6f %.6f] grados.\n', ...
    minimumCondition, alpha(minimumIndex), beta(minimumIndex), theta(minimumIndex));
fprintf('Error maximo de ortogonalidad de Q: %.3e | error de determinante: %.3e.\n', ...
    maximumOrthogonalityError, maximumDeterminantError);
fprintf('Error relativo de reconstruccion de covarianza de Cholesky: %.3e.\n', ...
    covarianceError);

function cellMap = occupied_condition_cells(alpha, beta, theta, values, limitDegrees, cellWidthDegrees)
    edges = -limitDegrees:cellWidthDegrees:limitDegrees;
    centers = edges(1:end-1) + cellWidthDegrees/2;
    binCount = numel(centers);
    alphaBin = discretize(alpha, edges);
    betaBin = discretize(beta, edges);
    thetaBin = discretize(theta, edges);
    valid = ~isnan(alphaBin) & ~isnan(betaBin) & ~isnan(thetaBin);
    linearBins = sub2ind([binCount, binCount, binCount], alphaBin(valid), ...
        betaBin(valid), thetaBin(valid));
    cellMinima = accumarray(linearBins, values(valid), ...
        [binCount^3, 1], @min, NaN);
    occupied = find(isfinite(cellMinima));
    [alphaIndex, betaIndex, thetaIndex] = ind2sub( ...
        [binCount, binCount, binCount], occupied);
    cellMap = struct('alpha', centers(alphaIndex).', ...
        'beta', centers(betaIndex).', 'theta', centers(thetaIndex).', ...
        'values', cellMinima(occupied));
end

function label = thousands_label(value)
    label = sprintf('%d', value);
    separatorPositions = numel(label)-2:-3:2;
    for position = separatorPositions
        label = [label(1:position-1), ',', label(position:end)]; %#ok<AGROW>
    end
end
