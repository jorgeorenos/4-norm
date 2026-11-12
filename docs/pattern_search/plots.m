% Reproduce the four Pattern Search winning-trajectory charts
% from the data saved in docs/patter_search/data.mat. Each chart plots
% the per-round estimate of kappa_4(PQ) for the trajectory that won among
% 100 Pattern Search starts, against the dotted Q = I reference
% kappa_4(P). 
scriptDirectory = fullfile('docs', 'patter_search');
dataFile = fullfile(scriptDirectory, 'data.mat');
if ~isfile(dataFile)
    error('plots:MissingData', ...
        'data.mat was not found next to plots.m: %s', dataFile);
end
data = load(dataFile);

charts = struct( ...
    'field', {'so2', 'so3', 'so5', 'so9'}, ...
    'fileName', {'pattern_search_cond4_2x2.png', ...
                 'pattern_search_cond4_3x3.png', ...
                 'pattern_search_cond4_5x5.png', ...
                 'pattern_search_cond4_9x9.png'});

for chartIndex = 1:numel(charts)
    fieldName = charts(chartIndex).field;
    dataset = data.(fieldName);
    plot_pattern_search_trajectory(dataset, ...
        fullfile(scriptDirectory, charts(chartIndex).fileName));
end

function plot_pattern_search_trajectory(dataset, outputFile)
% Draw one winning-trajectory chart from a harvested dataset struct and
% export it as a PNG at outputFile.
figureHandle = figure('Name', sprintf('Pattern Search winner in SO(%d)', ...
    dataset.K), 'Color', 'w', 'Position', [100, 100, 1000, 600]);
axesHandle = axes('Parent', figureHandle);
plot(axesHandle, dataset.outerIterations, dataset.historyValues, 'o-', ...
    'Color', [0.10, 0.25, 0.70], 'LineWidth', 1.4, ...
    'MarkerFaceColor', [0.10, 0.25, 0.70], ...
    'DisplayName', 'Winning trajectory among 100 starts');
hold(axesHandle, 'on');
yline(axesHandle, dataset.baselineValue, 'k:', 'LineWidth', 1.5, ...
    'DisplayName', 'Reference: kappa_4(P), Q = I');
hold(axesHandle, 'off');
grid(axesHandle, 'on');
xlabel(axesHandle, 'iteration');
ylabel(axesHandle, 'estimated kappa_4(PQ)');
title(axesHandle, sprintf(['SO(%d): winner %d (sample %d), ', ...
    '%.6f to %.6f'], dataset.K, dataset.bestStartIndex, ...
    dataset.selectedSampleIndex, dataset.initialValue, dataset.finalValue));
progressYLimits = ylim(axesHandle);
axisPadding = 0.05*diff(progressYLimits);
ylim(axesHandle, [progressYLimits(1) - axisPadding, ...
    progressYLimits(2) + axisPadding]);
legend(axesHandle, 'Location', 'best');

outputDirectory = fileparts(outputFile);
if ~isfolder(outputDirectory)
    mkdir(outputDirectory);
end
exportgraphics(figureHandle, outputFile, 'Resolution', 200);
close(figureHandle);
end
