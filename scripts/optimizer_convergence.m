%OPTIMIZER_CONVERGENCE Compara estimaciones multinicio de ||Q||_{4->4}.
%   Se puede ejecutar desde cualquier directorio. Para cada matriz Haar Q,
%   estima la norma inducida 4 con conjuntos anidados de inicios aleatorios.
%   Al restablecer el estado aleatorio antes de cada estimacion, se reutilizan
%   las primeras N direcciones de una misma secuencia, manteniendo Q fija.
%   La estimacion con mas inicios es una referencia numerica finita;
%   la brecha respecto de ella no certifica un maximo global.
%   Los resultados quedan en rotations, norm_q_values, convergence_values,
%   reference_values, relative_gaps y summary. Abre una figura de dos
%   paneles y requiere unicamente MATLAB base.

%% Localizar las funciones a partir de la ubicacion del script.
scriptDirectory = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDirectory);
addpath(fullfile(projectRoot, 'src', '4-norm'));
addpath(fullfile(projectRoot, 'src', 'helpers'));

%% Configuracion.
DIMENSION = 9;
N_MATRICES = 10;
START_COUNTS = [1, 2, 5, 10, 25, 50, 100, 250, 500, 750, 1000];
MAX_STARTS = START_COUNTS(end);
SEED = 7;
MAX_ITERATIONS = 300;
NORM_TOLERANCE = 1e-10;
REFERENCE_TOLERANCE = 1e-8;

rng(SEED, 'twister');
normHandle = @norm_4;
baseOptions = struct('MaxIterations', MAX_ITERATIONS, ...
    'NormTolerance', NORM_TOLERANCE);
numStartCounts = numel(START_COUNTS);

rotations = zeros(DIMENSION, DIMENSION, N_MATRICES);
norm_q_values = zeros(N_MATRICES, numStartCounts);

for matrix_index = 1:N_MATRICES
    Q = haar_so(DIMENSION);
    rotations(:, :, matrix_index) = Q;

    % Restablecer este estado conserva conjuntos anidados de inicios.
    startsState = rng;
    for start_count_index = 1:numStartCounts
        options = baseOptions;
        options.NumRandomStarts = START_COUNTS(start_count_index);

        rng(startsState);
        norm_q_values(matrix_index, start_count_index) = ...
            compute_induced_norm(Q, normHandle, options);
    end

    % Avanzar el estado como en una estimacion con MAX_STARTS direcciones.
    rng(startsState);
    generate_unit_l4_vectors(DIMENSION, MAX_STARTS, 'random');
end

% El mejor valor sobre conjuntos anidados no debe disminuir.
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

%% Mostrar la estimacion acumulada y la brecha frente a la referencia.
figure('Name', 'Convergencia multinicio de la norma inducida 4', 'Color', 'w');
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile
hold on
for matrix_index = 1:N_MATRICES
    semilogx(START_COUNTS, convergence_values(matrix_index, :), '-o', ...
        'DisplayName', sprintf('Q_%d', matrix_index));
end
set(gca, 'XScale', 'log')
grid on
xlabel('Numero de inicios aleatorios')
ylabel('Mejor estimacion de ||Q||_{4\rightarrow4}')
title('Estimacion multinicio de la norma inducida 4')
legend('Location', 'southeast')

nexttile
hold on
for matrix_index = 1:N_MATRICES
    % En escala logaritmica, representar las brechas nulas al nivel de eps.
    semilogy(START_COUNTS, max(relative_gaps(matrix_index, :), eps), '-o', ...
        'DisplayName', sprintf('Q_%d', matrix_index));
end
set(gca, 'YScale', 'log')
grid on
xlabel('Numero de inicios aleatorios')
ylabel(sprintf('Brecha relativa frente a %d inicios', MAX_STARTS))
title('Brecha frente a la referencia finita')
legend('Location', 'southwest')

sgtitle(sprintf('Convergencia de la norma inducida 4 para %d matrices Haar en SO(%d)', ...
    N_MATRICES, DIMENSION))
