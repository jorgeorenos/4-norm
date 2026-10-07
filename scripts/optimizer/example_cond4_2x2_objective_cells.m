% Visualize binned kappa_4(Q) estimates in the canonical SO(2) domain.
% Each point is one occupied angular cell. Its height and color are the
% minimum fixed-bank estimate found among the sampled rotations in that cell.
scriptDirectory = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptDirectory));
addpath(fullfile(projectRoot, 'src', '4-norm_handle'), '-begin');
addpath(fullfile(projectRoot, 'src', 'helpers'), '-begin');
addpath(fullfile(projectRoot, 'src', 'optimizer'), '-begin');

K = 2;
numSamples = 10000;
batchSize = 1000;
numRandomStarts = 75;
angleLimitDegrees = 45;
cellWidthDegrees = 0.25;
seeds = struct('Haar', 42, 'Starts', 43);

% P = I makes compute_cond4_pq_fixed_starts evaluate kappa_4(Q), while
% fixed forward/inverse banks make all sampled rotations comparable.
P = eye(K);
normOptions = struct('MaxIterations', 1000, ...
    'NormTolerance', 1e-10, 'StationarityTolerance', 1e-8, ...
    'FeasibilityTolerance', 1e-12, 'MaxWorkingMemoryMB', 128, ...
    'ScreenIterations', 2, 'NumFinalists', 4);
normHandle = @norm4_columns;
banks = create_norm4_start_banks(K, numRandomStarts, seeds.Starts);

rng(seeds.Haar, 'twister');
canonicalAngles = zeros(numSamples, 1);
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

    canonicalRotations = canonicalize_so2_signed_permutations(rotations);
    conditionNumbers(batchIndices) = compute_cond4_pq_fixed_starts( ...
        P, canonicalRotations, banks, normHandle, normOptions);
    canonicalAngles(batchIndices) = reshape(atan2(canonicalRotations(2,1,:), ...
        canonicalRotations(1,1,:)), [], 1);

    for page = 1:count
        Q = canonicalRotations(:,:,page);
        maximumOrthogonalityError = max(maximumOrthogonalityError, ...
            norm(Q'*Q - eye(K), 'fro'));
        maximumDeterminantError = max(maximumDeterminantError, ...
            abs(det(Q) - 1));
    end
end
elapsedSeconds = toc(timer);

assert(all(isfinite(conditionNumbers)) && all(conditionNumbers > 0));
assert(all(isfinite(canonicalAngles)));
assert(all(abs(canonicalAngles) <= pi/4 + 1e-12));
assert(maximumOrthogonalityError < 1e-12*max(1, K));
assert(maximumDeterminantError < 1e-12*max(1, K));

anglesDegrees = rad2deg(canonicalAngles);
[minimumCondition, minimumIndex] = min(conditionNumbers);
cellMap = occupied_condition_cells(anglesDegrees, conditionNumbers, ...
    angleLimitDegrees, cellWidthDegrees);
occupiedCount = numel(cellMap.values);

figureHandle = figure('Name', 'SO(2) celdas angulares ocupadas', 'Color', 'w');
axesHandle = axes(figureHandle);
% title
title('Evolucion del numero de condicion κ₄(Q)');
% subtitle
subtitle(sprintf('K=%d · %s rotaciones · Mejor \\kappa = %.6f', ...
    K, thousands_label(numSamples), minimumCondition));
hold(axesHandle, 'on');
scatter(axesHandle, cellMap.angle, cellMap.values, 22, log10(cellMap.values), ...
    'filled', 'MarkerEdgeColor', 'none', 'DisplayName', 'Numero de condicion κ₄(Q)');
plot(axesHandle, anglesDegrees(minimumIndex), minimumCondition, 'rp', ...
    'MarkerSize', 11, 'MarkerFaceColor', 'r', 'LineWidth', 1.0, ...
    'DisplayName', 'Minimo muestral');
hold(axesHandle, 'off');
colormap(figureHandle, turbo(256));
colorLimits = log10([min(cellMap.values), max(cellMap.values)]);
if diff(colorLimits) == 0
    colorLimits = colorLimits + [-0.01, 0.01];
end
clim(axesHandle, colorLimits);
colorbarHandle = colorbar(axesHandle);
ylabel(colorbarHandle, 'log_{10}(minimo estimado de \kappa_4(Q))');

grid(axesHandle, 'on');
xlim(axesHandle, [-angleLimitDegrees, angleLimitDegrees]);
ylim(axesHandle, [0.99, max([sqrt(2), max(cellMap.values)]) + 0.02]);
xlabel(axesHandle, '\theta (grados)');
ylabel(axesHandle, '\kappa_4(Q)');
legend(axesHandle, 'Location', 'southwest');

fprintf('\nMapa de celdas angulares ocupadas SO(2)\n');
fprintf('Se evaluaron %d matrices Q canonicas en %.3f segundos.\n', ...
    numSamples, elapsedSeconds);
fprintf('Se muestran %d de %d celdas angulares ocupadas.\n', ...
    occupiedCount, occupiedCount);
fprintf('Minimo muestral fijo de kappa_4(Q): %.12f en theta = %.6f grados.\n', ...
    minimumCondition, anglesDegrees(minimumIndex));
fprintf('Error maximo de ortogonalidad de Q: %.3e | error de determinante: %.3e.\n', ...
    maximumOrthogonalityError, maximumDeterminantError);

function cellMap = occupied_condition_cells(angles, values, limitDegrees, cellWidthDegrees)
    edges = -limitDegrees:cellWidthDegrees:limitDegrees;
    centers = edges(1:end-1) + cellWidthDegrees/2;
    angleBin = discretize(angles, edges);
    valid = ~isnan(angleBin);
    cellMinima = accumarray(angleBin(valid), values(valid), ...
        [numel(centers), 1], @min, NaN);
    occupied = find(isfinite(cellMinima));
    cellMap = struct('angle', centers(occupied).', 'values', cellMinima(occupied));
end

function label = thousands_label(value)
    label = sprintf('%d', value);
    separatorPositions = numel(label)-2:-3:2;
    for position = separatorPositions
        label = [label(1:position-1), ',', label(position:end)]; %#ok<AGROW>
    end
end
