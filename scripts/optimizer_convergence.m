%OPTIMIZER_CONVERGENCE Compare nested multistart estimates for four SO(9) rotations.
%   Run this script from any current directory. It changes to the OSST
%   repository root, activates the project, draws four reproducible Haar
%   rotations in SO(9), and estimates each
%   induced 4-norm with nested sets of random initial vectors. The estimate
%   at N starts always reuses the first N vectors from the same pre-generated
%   sequence, so its best-so-far value cannot decrease as N increases.
%
%   For each rotation Q, the nonlinear power iteration approximates
%
%       ||Q||_{4->4} = max_{||x||_4=1} ||Q*x||_4.
%
%   The final figure has two panels: the best norm estimate and its relative
%   gap to the estimate based on MAX_STARTS starts. The latter is a numerical
%   reference, not a certificate of the global maximum. Results remain in
%   the workspace as rotations, best_values, convergence_values, and summary.
%   The script advances the global RNG, prints a summary table, and opens one
%   figure; it does not write files. Requires only base MATLAB and OSST.
%
%   See also ORTHOGONAL_MATRIX_GENERATOR, RNG, NORM.

%% Configuration.
DIMENSION = 30;
NORM_ORDER = 4;
N_MATRICES = 6;
START_COUNTS = [1, 2, 5, 10, 25, 50, 100, 250, 500, 750, 1000];
MAX_STARTS = START_COUNTS(end);
SEED = 7;
MAX_ITERATIONS = 300;
RELATIVE_TOLERANCE = 1e-10;
REFERENCE_TOLERANCE = 1e-8;

% One stream generates both the rotations and their nested initial vectors.
rng(SEED, 'twister');
rotations = zeros(DIMENSION, DIMENSION, N_MATRICES);
best_values = zeros(N_MATRICES, MAX_STARTS);

for matrix_index = 1:N_MATRICES
    [Q, upper_triangular] = qr(randn(DIMENSION), 0);
    diagonal_signs = sign(diag(upper_triangular));
    diagonal_signs(diagonal_signs == 0) = 1;
    Q = Q*diag(diagonal_signs);
    if det(Q) < 0
        Q(:, end) = -Q(:, end);
    end
    rotations(:, :, matrix_index) = Q;

    % Columns are fixed before optimization: start-count comparisons are fair.
    initial_points = randn(DIMENSION, MAX_STARTS);
    initial_norms = sum(abs(initial_points).^NORM_ORDER, 1).^(1/NORM_ORDER);
    initial_points = bsxfun(@rdivide, initial_points, initial_norms);
    best_value = -Inf;
    for start_index = 1:MAX_STARTS
        x = initial_points(:, start_index);
        candidate_value = norm(Q*x, NORM_ORDER);
        previous_value = candidate_value;
        dual_exponent_minus_one = 1/(NORM_ORDER - 1);
        for iteration_index = 1:MAX_ITERATIONS
            y = Q*x;
            gradient = Q'*(abs(y).^(NORM_ORDER - 2).*y);
            x = sign(gradient).*abs(gradient).^dual_exponent_minus_one;
            x = x/norm(x, NORM_ORDER);

            current_value = norm(Q*x, NORM_ORDER);
            candidate_value = max(candidate_value, current_value);
            if abs(current_value - previous_value) <= ...
                    RELATIVE_TOLERANCE*max(1, current_value)
                break
            end
            previous_value = current_value;
        end
        best_value = max(best_value, candidate_value);
        best_values(matrix_index, start_index) = best_value;
    end
end

convergence_values = best_values(:, START_COUNTS);
reference_values = best_values(:, MAX_STARTS);
relative_gaps = (reference_values - convergence_values) ./ reference_values;
converged_start_counts = zeros(N_MATRICES, 1);
for matrix_index = 1:N_MATRICES
    first_converged_index = find( ...
        relative_gaps(matrix_index, :) <= REFERENCE_TOLERANCE, 1, 'first');
    converged_start_counts(matrix_index) = START_COUNTS(first_converged_index);
end

matrix_id = (1:N_MATRICES).';
summary = table(matrix_id, convergence_values(:, 6), reference_values, ...
    relative_gaps(:, 6), converged_start_counts, ...
    'VariableNames', {'matrix', 'estimate_50_starts', 'reference_estimate', ...
    'relative_gap_50_starts', 'first_reference_tolerance_start_count'});
disp(summary)

%% Plot the nested multistart convergence for every rotation.
figure('Name', 'Induced 4-norm multistart convergence', 'Color', 'w');
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile
hold on
for matrix_index = 1:N_MATRICES
    semilogx(START_COUNTS, convergence_values(matrix_index, :), '-o', ...
        'DisplayName', sprintf('Q_%d', matrix_index));
end
grid on
xlabel('Numero de inicios aleatorios')
ylabel('Mejor estimación de la norma inducida 4')
title('Mejor estimación de la norma 4 inducida')
legend('Location', 'southeast')

nexttile
hold on
for matrix_index = 1:N_MATRICES
    semilogy(START_COUNTS, max(relative_gaps(matrix_index, :), eps), '-o', ...
        'DisplayName', sprintf('Q_%d', matrix_index));
end
grid on
xlabel('Numero de inicios aleatorios')
ylabel('Brecha relativa con respecto a la referencia de 1000 inicios')
title('Diagnóstico de convergencia')
legend('Location', 'southwest')

sgtitle(sprintf(['Convergencia de la norma inducida %d para %d Matrices Q en SO(%d) '], ...
    NORM_ORDER, N_MATRICES, DIMENSION))
