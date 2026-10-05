%RENDER_CONVERGENCE_ANALYSIS Render convergence figures for docs/optimizer.qmd.
%   Run from any current directory. The script writes one two-panel figure
%   per dimension to docs/optimizer/figures/convergence_analysis. Each curve estimates
%   the induced 4-norm ||Q||_{4->4} using nested sets of random starts.

scriptDirectory = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptDirectory));
outputDir = fullfile(scriptDirectory, 'figures', 'convergence_analysis');
if ~isfolder(outputDir)
    mkdir(outputDir);
end

addpath(fullfile(projectRoot, 'src', '4-norm'));
addpath(fullfile(projectRoot, 'src', 'helpers'));

dimensions = [2, 5, 9, 17];
numMatrices = 10;
startCounts = [1, 2, 5, 10, 25, 50, 100, 250, 500, 750, 1000];
maxIterations = 1000;
normTolerance = 1e-10;
seed = 7;

for dimensionIndex = 1:numel(dimensions)
    dimension = dimensions(dimensionIndex);
    [inducedNormValues, relativeGaps] = ...
        estimate_convergence(dimension, numMatrices, startCounts, ...
        maxIterations, normTolerance, seed);

    fig = plot_convergence(dimension, startCounts, inducedNormValues, ...
        relativeGaps);
    filename = sprintf('convergence-analysis-so%d.png', dimension);
    % Export the figure to the output directory
    exportgraphics(fig, fullfile(outputDir, filename), 'Resolution', 150);
    close(fig);
end

function [inducedNormValues, relativeGaps] = ...
        estimate_convergence(dimension, numMatrices, startCounts, ...
        maxIterations, normTolerance, seed)
    rng(seed, 'twister');
    numStartCounts = numel(startCounts);
    inducedNormValues = zeros(numMatrices, numStartCounts);
    normHandle = @norm_4;
    baseOptions = struct('MaxIterations', maxIterations, ...
        'NormTolerance', normTolerance);

    for matrixIndex = 1:numMatrices
        Q = haar_so(dimension);
        startsState = rng;
        for startCountIndex = 1:numStartCounts
            options = baseOptions;
            options.NumRandomStarts = startCounts(startCountIndex);

            rng(startsState);
            inducedNormValues(matrixIndex, startCountIndex) = ...
                compute_induced_norm(Q, normHandle, options);
        end

        % Keep the rotation stream aligned with the maximum-start run.
        rng(startsState);
        directions = randn(dimension, startCounts(end));
        zeroColumns = ~any(directions, 1);
        while any(zeroColumns)
            directions(:,zeroColumns) = randn(dimension, sum(zeroColumns));
            zeroColumns = ~any(directions, 1);
        end
    end

    assert(all(all(diff(inducedNormValues, 1, 2) >= -128*eps)), ...
        'Nested multistart induced-norm estimates must not decrease.');
    referenceValues = inducedNormValues(:, end);
    relativeGaps = (referenceValues - inducedNormValues) ./ referenceValues;
end

function fig = plot_convergence(dimension, startCounts, inducedNormValues, ...
        relativeGaps)
    fig = figure('Visible', 'off', 'Color', 'w', ...
        'Position', [100, 100, 1200, 460]);
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile
    hold on
    for matrixIndex = 1:size(inducedNormValues, 1)
        semilogx(startCounts, inducedNormValues(matrixIndex, :), '-o', ...
            'DisplayName', sprintf('Q_%d', matrixIndex));
    end
    set(gca, 'XScale', 'log')
    grid on
    xlabel('Numero de inicios aleatorios')
    ylabel('Mejor estimacion de ||Q||_{4->4}')
    title('Estimacion multinicio')
    legend('Location', 'best')

    nexttile
    hold on
    for matrixIndex = 1:size(relativeGaps, 1)
        semilogy(startCounts, max(relativeGaps(matrixIndex, :), eps), '-o', ...
            'DisplayName', sprintf('Q_%d', matrixIndex));
    end
    grid on
    xlabel('Numero de inicios aleatorios')
    ylabel('Brecha relativa frente a 1000 inicios')
    title('Brecha frente a la referencia finita')
    legend('Location', 'best')

    sgtitle(sprintf(['Convergencia de la norma inducida 4 para 10 ' ...
        'matrices Haar en SO(%d)'], ...
        dimension))
end
