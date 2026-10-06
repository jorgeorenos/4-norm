%OPTIMIZER_CONVERGENCE Compares multistart estimates of ||Q||_{4->4}.
%   Can be run from any directory. For each Haar matrix Q, estimates the
%   induced 4-norm using nested sets of random starting points. Resetting the
%   random state before each estimate reuses the first N directions from the
%   same sequence while keeping Q fixed. The estimate with the most starts is
%   a finite numerical reference; the gap to it does not certify a global
%   maximum. Results are stored in rotations, norm_q_values,
%   convergence_values, reference_values, relative_gaps, and summary. Opens a
%   two-panel figure and requires MATLAB base only.

%% Locate functions relative to the script location.
scriptDirectory = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDirectory);
addpath(fullfile(projectRoot, 'src', '4-norm'));
addpath(fullfile(projectRoot, 'src', 'helpers'));

%% Configuration.
DIMENSION = 9;
N_MATRICES = 10;
START_COUNTS = [1, 2, 5, 10, 25, 50, 100, 250, 500, 750, 1000];
MAX_STARTS = START_COUNTS(end);
SEED = 7;
MAX_ITERATIONS = 300;
NORM_TOLERANCE = 1e-10;
REFERENCE_TOLERANCE = 1e-8;

rng(SEED, 'twister');
baseOptions = struct('MaxIterations', MAX_ITERATIONS, ...
    'NormTolerance', NORM_TOLERANCE, ...
    'StationarityTolerance', 1e-8, 'FeasibilityTolerance', 1e-12);
numStartCounts = numel(START_COUNTS);

rotations = zeros(DIMENSION, DIMENSION, N_MATRICES);
norm_q_values = zeros(N_MATRICES, numStartCounts);

for matrix_index = 1:N_MATRICES
    Q = haar_so(DIMENSION);
    rotations(:, :, matrix_index) = Q;

    % A single batched iteration preserves the same starting directions for
    % every start count. Best prefixes replace 11 separate calls without
    % changing matrix Q or the state of the random sequence.
    starts = randn(DIMENSION, MAX_STARTS);
    zeroColumns = ~any(starts, 1);
    while any(zeroColumns)
        starts(:,zeroColumns) = randn(DIMENSION, sum(zeroColumns));
        zeroColumns = ~any(starts, 1);
    end
    perStartValues = power_norm4_multiple_starts(Q, starts, baseOptions);
    bestPrefix = cummax(perStartValues(:));
    norm_q_values(matrix_index, :) = bestPrefix(START_COUNTS).';
end

% The best value over nested sets must not decrease.
assert(all(all(diff(norm_q_values, 1, 2) >= -128*eps)), ...
    'Nested multistart induced-norm estimates must not decrease.');

convergence_values = norm_q_values;
reference_values = norm_q_values(:, end);
relative_gaps = (reference_values - convergence_values) ./ reference_values;
converged_start_counts = zeros(N_MATRICES, 1);
for matrix_index = 1:N_MATRICES
    first_converged_index = find( ...
        relative_gaps(matrix_index, :) <= REFERENCE_TOLERANCE, 1, 'first');
    converged_start_counts(matrix_index) = START_COUNTS(first_converged_index);
end

reportStartCount = 50;
reportStartIndex = find(START_COUNTS == reportStartCount, 1, 'first');
matrix_id = (1:N_MATRICES).';
summary = table(matrix_id, norm_q_values(:, reportStartIndex), ...
    reference_values, ...
    relative_gaps(:, reportStartIndex), converged_start_counts, ...
    'VariableNames', {'matrix', 'norm_q_50_starts', ...
    'reference_norm_q', 'relative_gap_50_starts', ...
    'first_reference_tolerance_start_count'});
disp(summary)

%% Plot the cumulative estimate and the gap to the reference.
figure('Name', 'Multistart convergence of the induced 4-norm', 'Color', 'w');
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile
hold on
for matrix_index = 1:N_MATRICES
    semilogx(START_COUNTS, convergence_values(matrix_index, :), '-o', ...
        'DisplayName', sprintf('Q_%d', matrix_index));
end
set(gca, 'XScale', 'log')
grid on
xlabel('Number of random starts')
ylabel('Best estimate of ||Q||_{4\rightarrow4}')
title('Multistart estimate of the induced 4-norm')
legend('Location', 'southeast')

nexttile
hold on
for matrix_index = 1:N_MATRICES
    % On a logarithmic scale, plot zero gaps at the eps level.
    semilogy(START_COUNTS, max(relative_gaps(matrix_index, :), eps), '-o', ...
        'DisplayName', sprintf('Q_%d', matrix_index));
end
set(gca, 'YScale', 'log')
grid on
xlabel('Number of random starts')
ylabel(sprintf('Relative gap to %d starts', MAX_STARTS))
title('Gap to the finite reference')
legend('Location', 'southwest')

sgtitle(sprintf('Induced 4-norm convergence for %d Haar matrices in SO(%d)', ...
    N_MATRICES, DIMENSION))
