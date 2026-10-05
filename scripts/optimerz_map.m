% Compare scalar power estimates from fixed starts at increasing budgets.
scriptDirectory = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDirectory);
addpath(fullfile(projectRoot, 'src', '4-norm'));
addpath(fullfile(projectRoot, 'src', 'helpers'));
rng(42, 'twister');
Q = haar_so(2);
starts = randn(2, 10);
X0 = starts ./ norm_4(starts, 1);
iterationBudgets = [1 2 5 10 20 50 100 250];
options = struct('NormTolerance', 1e-10, ...
    'StationarityTolerance', 1e-8, 'FeasibilityTolerance', 1e-12);
values = zeros(size(X0,2), numel(iterationBudgets));
for k = 1:numel(iterationBudgets)
    options.MaxIterations = iterationBudgets(k);
    % Keep all ten fixed starts in one vectorized power iteration.
    values(:,k) = power_norm4_multiple_starts(Q, @norm_4, X0, options);
end
figure;
semilogx(iterationBudgets, values.', '-o');
grid on;
xlabel('Limite de iteraciones');
ylabel('Mejor amplitud factible ||Qx||_4');
title('Estimaciones desde diez inicios fijos');
% Each point restarts from the same x0; it is not a stored iteration history.
