% Reproduce the four optimizer geometry/result charts in English from
% docs/optimizer/data_geometry_analysis.mat. This script only redraws
% saved data; see generate_geometry_data.m for how that data was
% computed (mirrors render_optimizer_figures.m).
scriptDirectory = fileparts(mfilename('fullpath'));
dataFile = fullfile(scriptDirectory, 'data_geometry_analysis.mat');
if ~isfile(dataFile)
    error('plot_geometry_analysis:MissingData', ...
        'data_geometry_analysis.mat was not found next to this script: %s', dataFile);
end
data = load(dataFile);

outputDir = fullfile(scriptDirectory, 'figures', 'generated');
if ~isfolder(outputDir)
    mkdir(outputDir);
end

% --- S_4 geometry and amplitude in SO(2) --------------------------------
d = data.geometry2d;
Xclosed = [d.X, d.X(:,1)];
Yclosed = [d.Y, d.Y(:,1)];
allCoords = [Xclosed, Yclosed];
coordLimit = max(abs(allCoords(:)));
margin = 0.08 * coordLimit;
plotLimits = [-coordLimit - margin, coordLimit + margin];

fig = figure('Visible', 'off', 'Position', [100 100 1400 400]);
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot(Xclosed(1,:), Xclosed(2,:), 'LineWidth', 1.8);
axis equal;
xlim(plotLimits); ylim(plotLimits);
grid on;
xlabel('x_1'); ylabel('x_2');
title('S_4 unit sphere');
nexttile;
plot(Xclosed(1,:), Xclosed(2,:), '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.5);
hold on;
plot(Yclosed(1,:), Yclosed(2,:), 'Color', [0.85 0.33 0.10], 'LineWidth', 1.8);
axis equal;
xlim(plotLimits); ylim(plotLimits);
grid on;
xlabel('u_1'); ylabel('u_2');
title('S_4 and its image Q(S_4)');
legend('Original: ||x||_4=1', 'Transformed: y=Qx', 'Location', 'best');
hold off;
nexttile;
plot(d.theta, d.amplitudes, 'Color', [0.15 0.35 0.70], 'LineWidth', 1.2);
grid on;
xlabel('\theta (radians)');
ylabel('||Qx(\theta)||_4');
title('Amplitude over Q(S_4)');
xlim([0 2*pi]);
exportgraphics(fig, fullfile(outputDir, 'optimizer-geometry-2d.png'), 'Resolution', 150);
close(fig);

% --- S_4 surface colored by amplitude in SO(3) --------------------------
d = data.geometry3d;
X1 = reshape(d.Xgrid(1,:), size(d.azGrid));
X2 = reshape(d.Xgrid(2,:), size(d.azGrid));
X3 = reshape(d.Xgrid(3,:), size(d.azGrid));
fig = figure('Visible', 'off', 'Position', [100 100 800 600]);
surf(X1, X2, X3, d.amplitudes, 'EdgeColor', 'none', 'FaceColor', 'interp');
axis equal; grid on;
xlabel('x_1'); ylabel('x_2'); zlabel('x_3');
title('S_4 colored by ||Qx||_4');
colorbar;
view(3);
exportgraphics(fig, fullfile(outputDir, 'optimizer-geometry-3d.png'), 'Resolution', 150);
close(fig);

% --- Induced 4-norm estimate against an angular scan in SO(2) -----------
d = data.result2d;
fig = figure('Visible', 'off', 'Position', [100 100 1000 500]);
plot([d.theta, 2*pi], [d.amplitudes, d.amplitudes(1)], ...
    'Color', [0.15 0.35 0.70], 'LineWidth', 1.5);
hold on;
yline(d.qNorm4, '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.5);
xlim([0 2*pi]); grid on;
xlabel('\theta (radians)');
ylabel('||Qx(\theta)||_4');
title(sprintf('Q = [%.3f, %.3f; %.3f, %.3f] | ||Q||_4 \\approx %.6f', ...
    d.Q(1,1), d.Q(1,2), d.Q(2,1), d.Q(2,2), d.qNorm4));
legend('Angular grid', 'Induced 4-norm estimate', 'Location', 'best');
hold off;
exportgraphics(fig, fullfile(outputDir, 'optimizer-result-2d.png'), 'Resolution', 150);
close(fig);

% --- Induced 4-norm estimate against a spherical grid in SO(3) ----------
d = data.result3d;
X1 = reshape(d.Xgrid(1,:), size(d.azGrid));
X2 = reshape(d.Xgrid(2,:), size(d.azGrid));
X3 = reshape(d.Xgrid(3,:), size(d.azGrid));
fig = figure('Visible', 'off', 'Position', [100 100 900 700]);
surf(X1, X2, X3, d.amplitudes, 'EdgeColor', 'none', 'FaceColor', 'interp');
hold on;
plot3(d.xRef(1), d.xRef(2), d.xRef(3), 'p', 'Color', [0.85 0.33 0.10], ...
    'MarkerFaceColor', [0.85 0.33 0.10], 'MarkerSize', 12);
axis equal; grid on;
xlabel('x_1'); ylabel('x_2'); zlabel('x_3');
title(sprintf('S_4 colored by ||Qx||_4 (n=3) | ||Q||_4 \\approx %.6f', d.qNorm4));
legend('Surface', 'Refined maximum (grid)', 'Location', 'best');
colorbar; view(3);
hold off;
exportgraphics(fig, fullfile(outputDir, 'optimizer-result-3d.png'), 'Resolution', 150);
close(fig);

fprintf('Wrote 4 geometry/result charts to %s\n', outputDir);
