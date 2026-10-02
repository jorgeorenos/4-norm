% Geometric example of the 4-norm unit sphere in R^2.
% A rotation in SO(2) preserves Euclidean lengths, angles, and area, but
% generally not the 4-norm; therefore Y is not renormalized here.

scriptPath = mfilename('fullpath');
scriptDirectory = fileparts(scriptPath);
projectRoot = fileparts(scriptDirectory);
addpath(fullfile(projectRoot, 'src', '4-norm'));
addpath(fullfile(projectRoot, 'src', 'helpers'));

rng(42, 'twister');
n = 2;
N = 1000;

X = generate_unit_l4_vectors(n, N, 'angular');
A = haar_so(n);
disp('Matrix A sampled from SO(2):');
disp(A);
Y = A * X;

xNorms = sum(abs(X).^4, 1).^(1/4);
xNormError = max(abs(xNorms - 1));
orthogonalityError = norm(A' * A - eye(n), 'fro');
determinantError = abs(det(A) - 1);
fprintf('Maximum 4-norm error for X: %.3e\n', xNormError);
fprintf('Orthogonality error of A: %.3e\n', orthogonalityError);
fprintf('Error of det(A) relative to 1: %.3e\n', determinantError);

Xclosed = [X, X(:, 1)];
Yclosed = [Y, Y(:, 1)];
allCoordinates = [Xclosed, Yclosed];
coordinateLimit = max(abs(allCoordinates(:)));
margin = 0.08 * coordinateLimit;
plotLimits = [-coordinateLimit - margin, coordinateLimit + margin];

figure;
tiledlayout(1, 2);

nexttile;
plot(Xclosed(1, :), Xclosed(2, :), 'LineWidth', 1.8);
axis equal;
xlim(plotLimits);
ylim(plotLimits);
grid on;
xlabel('x_1');
ylabel('x_2');
title('Unit sphere of the 4-norm');

nexttile;
plot(Xclosed(1, :), Xclosed(2, :), '--', 'Color', [0.45 0.45 0.45], ...
    'LineWidth', 1.5);
hold on;
plot(Yclosed(1, :), Yclosed(2, :), 'Color', [0.85 0.33 0.10], ...
    'LineWidth', 1.8);
axis equal;
xlim(plotLimits);
ylim(plotLimits);
grid on;
xlabel('u_1');
ylabel('u_2');
title('Unit sphere and its image under A');
legend('Original: 4-norm = 1', 'Transformed: y = Ax', ...
    'Location', 'best');
hold off;

% The figure illustrates geometry; it neither identifies x* nor computes
% an induced matrix norm. Rotations by multiples of pi/2 do preserve S4,
% but Haar sampling is not repeated to force that special case.
