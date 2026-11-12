% Reproduce the four convergence-analysis charts in English from
% docs/optimizer/data_convergence_analysis.mat. This script only redraws
% saved data; see generate_convergence_data.m for how that data was
% computed (mirrors render_convergence_analysis.m). The original charts
% use Spanish axis labels and titles; this script redraws the same
% curves with English text.
scriptDirectory = fileparts(mfilename('fullpath'));
dataFile = fullfile(scriptDirectory, 'data_convergence_analysis.mat');
if ~isfile(dataFile)
    error('plot_convergence_analysis:MissingData', ...
        'data_convergence_analysis.mat was not found next to this script: %s', dataFile);
end
data = load(dataFile);

outputDir = fullfile(scriptDirectory, 'figures', 'convergence_analysis');
if ~isfolder(outputDir)
    mkdir(outputDir);
end

fields = {'so2', 'so5', 'so9', 'so17'};
for fieldIndex = 1:numel(fields)
    dataset = data.(fields{fieldIndex});
    fig = plot_convergence(dataset.dimension, dataset.startCounts, ...
        dataset.inducedNormValues, dataset.relativeGaps);
    filename = sprintf('convergence-analysis-so%d.png', dataset.dimension);
    exportgraphics(fig, fullfile(outputDir, filename), 'Resolution', 150);
    close(fig);
end
fprintf('Wrote 4 convergence charts to %s\n', outputDir);

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
    xlabel('Number of random starts')
    ylabel('Best estimate of ||Q||_{4->4}')
    title('Multistart estimate')
    legend('Location', 'best')

    nexttile
    hold on
    for matrixIndex = 1:size(relativeGaps, 1)
        semilogy(startCounts, max(relativeGaps(matrixIndex, :), eps), '-o', ...
            'DisplayName', sprintf('Q_%d', matrixIndex));
    end
    grid on
    xlabel('Number of random starts')
    ylabel('Relative gap versus 1000 starts')
    title('Gap versus the finite reference')
    legend('Location', 'best')

    sgtitle(sprintf(['Convergence of the induced 4-norm for 10 ' ...
        'Haar matrices in SO(%d)'], ...
        dimension))
end
