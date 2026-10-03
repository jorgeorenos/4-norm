% Render the static figures used by docs/optimizer/optimizer.qmd.
scriptDirectory = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptDirectory));

outputDir = fullfile(projectRoot, 'docs', 'optimizer', 'figures', 'generated');
if ~isfolder(outputDir)
    mkdir(outputDir);
end

addpath(fullfile(projectRoot, 'src', '4-norm'));
addpath(fullfile(projectRoot, 'src', 'helpers'));

rng(42, 'twister');
n = 2;
N = 1000;
Q = haar_so(n);
X = generate_unit_l4_vectors(n, N, 'angular');
Y = Q * X;
Xclosed = [X, X(:,1)];
Yclosed = [Y, Y(:,1)];
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
theta = (0:N-1) * (2*pi/N);
amplitudes = zeros(1,N);
for k = 1:N
    amplitudes(k) = norm_4(Q * X(:,k));
end
nexttile;
plot(theta, amplitudes, 'Color', [0.15 0.35 0.70], 'LineWidth', 1.2);
grid on;
xlabel('\theta (radians)');
ylabel('||Qx(\theta)||_4');
title('Amplitude over Q(S_4)');
xlim([0 2*pi]);
exportgraphics(fig, fullfile(outputDir, 'optimizer-geometry-2d.png'), 'Resolution', 150);
close(fig);

rng(42, 'twister');
n = 3;
Q = haar_so(n);
[azGrid, ~, Xgrid, amplitudes] = sphericalGrid(Q);
X1 = reshape(Xgrid(1,:), size(azGrid));
X2 = reshape(Xgrid(2,:), size(azGrid));
X3 = reshape(Xgrid(3,:), size(azGrid));
fig = figure('Visible', 'off', 'Position', [100 100 800 600]);
surf(X1, X2, X3, amplitudes, 'EdgeColor', 'none', 'FaceColor', 'interp');
axis equal; grid on;
xlabel('x_1'); ylabel('x_2'); zlabel('x_3');
title('S_4 colored by ||Qx||_4');
colorbar;
view(3);
exportgraphics(fig, fullfile(outputDir, 'optimizer-geometry-3d.png'), 'Resolution', 150);
close(fig);

rng(42, 'twister');
n = 2;
Q = haar_so(n);
normHandle = @norm_4;
options = struct('NumRandomStarts', 50, 'MaxIterations', 1000);
qNorm4 = compute_induced_norm(Q, normHandle, options);
N = 10000;
theta = (0:N-1) * (2*pi/N);
X = generate_unit_l4_vectors(2, N, 'angular');
amplitudes = zeros(1,N);
for k = 1:N
    amplitudes(k) = normHandle(Q * X(:,k));
end
fig = figure('Visible', 'off', 'Position', [100 100 1000 500]);
plot([theta, 2*pi], [amplitudes, amplitudes(1)], 'Color', [0.15 0.35 0.70], 'LineWidth', 1.5);
hold on;
yline(qNorm4, '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.5);
xlim([0 2*pi]); grid on;
xlabel('\theta (radians)');
ylabel('||Qx(\theta)||_4');
title(sprintf('Q = [%.3f, %.3f; %.3f, %.3f] | ||Q||_4 \\approx %.6f', Q(1,1), Q(1,2), Q(2,1), Q(2,2), qNorm4));
legend('Angular grid', 'Induced 4-norm estimate', 'Location', 'best');
hold off;
exportgraphics(fig, fullfile(outputDir, 'optimizer-result-2d.png'), 'Resolution', 150);
close(fig);

rng(42, 'twister');
n = 3;
Q = haar_so(n);
normHandle = @norm_4;
options = struct('NumRandomStarts', 100, 'MaxIterations', 1000);
qNorm4 = compute_induced_norm(Q, normHandle, options);
[azGrid, polGrid, Xgrid, amplitudes] = sphericalGrid(Q);
[~, sortedIdx] = sort(amplitudes(:), 'descend');
refinedOptions = optimset('Display', 'off', 'MaxFunEvals', 1000, 'MaxIter', 1000, 'TolX', 1e-12, 'TolFun', 1e-12);
bestRefined = -Inf;
bestAngles = [];
for j = 1:min(12, numel(sortedIdx))
    idx = sortedIdx(j);
    init = [azGrid(idx); polGrid(idx)];
    [ang, negVal] = fminsearch(@(a) -spherical_amplitude(a, Q, normHandle), init, refinedOptions);
    if -negVal > bestRefined
        bestRefined = -negVal;
        bestAngles = ang;
    end
end
refDir = [cos(bestAngles(1))*sin(bestAngles(2)); sin(bestAngles(1))*sin(bestAngles(2)); cos(bestAngles(2))];
xRef = refDir / normHandle(refDir);
X1 = reshape(Xgrid(1,:), size(azGrid));
X2 = reshape(Xgrid(2,:), size(azGrid));
X3 = reshape(Xgrid(3,:), size(azGrid));
fig = figure('Visible', 'off', 'Position', [100 100 900 700]);
surf(X1, X2, X3, amplitudes, 'EdgeColor', 'none', 'FaceColor', 'interp');
hold on;
plot3(xRef(1), xRef(2), xRef(3), 'p', 'Color', [0.85 0.33 0.10], 'MarkerFaceColor', [0.85 0.33 0.10], 'MarkerSize', 12);
axis equal; grid on;
xlabel('x_1'); ylabel('x_2'); zlabel('x_3');
title(sprintf('S_4 colored by ||Qx||_4 (n=3) | ||Q||_4 \\approx %.6f', qNorm4));
legend('Surface', 'Refined maximum (grid)', 'Location', 'best');
colorbar; view(3);
hold off;
exportgraphics(fig, fullfile(outputDir, 'optimizer-result-3d.png'), 'Resolution', 150);
close(fig);
