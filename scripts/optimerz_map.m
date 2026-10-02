% optimerz_map.m
% Map the generalized power optimizer iteration for a matrix Q in SO(2).
% Display the objective function ||Qx||_4 at every iteration from several
% initial points on the 4-norm unit sphere.
%
% Generate two figures:
%   1. Trajectory from one x0.
%   2. A 2-by-5 tiled layout from 10 different x0 values.

scriptPath = mfilename('fullpath');
scriptDir = fileparts(scriptPath);
projectRoot = fileparts(scriptDir);
addpath(fullfile(projectRoot, 'src', '4-norm'));
addpath(fullfile(projectRoot, 'src', 'helpers'));

rng(42, 'twister');

n = 2;
Q = haar_so(n);
normHandle = @norm_4;

options = struct('MaxIterations', 1000, ...
    'NormTolerance', 1e-10, ...
    'StationarityTolerance', 1e-8, ...
    'FeasibilityTolerance', 1e-12, ...
    'StoreHistory', true);

fprintf('Matrix Q in SO(2):\n');
disp(Q);

%% Figure 1: evolution from one initial point
x0_single = generate_unit_l4_vectors(n, 1, 'random');
[value_single, xBest_single, runInfo_single] = power_norm4_single_start(Q, normHandle, x0_single, options);

figure('Name', 'Evolution from one x0');
iter_single = 0:(length(runInfo_single.history) - 1);
plot(iter_single, runInfo_single.history, '-o', 'LineWidth', 1.5, 'MarkerSize', 5, ...
    'Color', [0.15 0.35 0.70]);
hold on;
plot(runInfo_single.iterations, value_single, 'r*', 'MarkerSize', 12, 'LineWidth', 2);
hold off;
grid on;
xlabel('Iteration number');
ylabel('||Qx||_4');
title(sprintf(['x0 = [%.4f; %.4f]  |  x* = [%.4f; %.4f]\n' ...
    'Converged: %d  |  Iterations: %d  |  Final value: %.6f'], ...
    x0_single(1), x0_single(2), xBest_single(1), xBest_single(2), ...
    runInfo_single.converged, runInfo_single.iterations, value_single));
legend('Objective trajectory', 'Final value', 'Location', 'best');

%% Figure 2: evolution from 10 initial points (2-by-5)
N_starts = 10;
X0 = generate_unit_l4_vectors(n, N_starts, 'random');

figure('Name', 'Evolution from 10 x0 values');
tiledlayout(2, 5, 'TileSpacing', 'compact', 'Padding', 'compact');

for j = 1:N_starts
    x0 = X0(:, j);
    [value_j, xBest_j, runInfo_j] = power_norm4_single_start(Q, normHandle, x0, options);

    nexttile;
    iter_j = 0:(length(runInfo_j.history) - 1);
    plot(iter_j, runInfo_j.history, '-o', 'LineWidth', 1.2, 'MarkerSize', 3, ...
        'Color', [0.15 0.35 0.70]);
    hold on;
    plot(runInfo_j.iterations, value_j, 'r*', 'MarkerSize', 10, 'LineWidth', 1.5);
    hold off;
    grid on;
    xlabel('Iteration');
    ylabel('||Qx||_4');
    title(sprintf('x0=[%.3f;%.3f] | x*=[%.3f;%.3f] | it=%d', ...
        x0(1), x0(2), xBest_j(1), xBest_j(2), runInfo_j.iterations), 'FontSize', 8);
end

sgtitle(sprintf('Optimizer evolution from 10 initial points\nQ = [%.4f, %.4f; %.4f, %.4f]', ...
    Q(1,1), Q(1,2), Q(2,1), Q(2,2)));
