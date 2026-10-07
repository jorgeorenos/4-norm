% Visualize low condition-number representatives in the canonical SO(3) domain.
% The plot uses principal rotation-vector coordinates log(Q) and colors the
% selected representatives by their raw, pre-canonical condition estimates.
scriptDirectory = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptDirectory));
addpath(fullfile(projectRoot, 'src', '4-norm_handle'), '-begin');
addpath(fullfile(projectRoot, 'src', 'helpers'), '-begin');

n = 3;
numSamples = 50000;
numRetained = 1000;
batchSize = 100;
numRandomStarts = 75;
options = struct('NumRandomStarts', numRandomStarts, 'MaxIterations', 1000, ...
    'MaxWorkingMemoryMB', 128, 'ScreenIterations', 2, 'NumFinalists', 4);
normHandle = @norm4_columns;

rng(42, 'twister');
conditionNumberFn = @(rotations) evaluate_condition_number_batch( ...
    rotations, normHandle, options);
[canonicalRotations, selectedConditionNumbers, selectedIndices] = ...
    initialize_condition_number(n, numSamples, batchSize, ...
    conditionNumberFn, numRetained);

% For Q = exp([omega]_x), the principal logarithm is [omega]_x. The
% resulting rotation vector omega is a local three-coordinate chart around
% the identity-centered canonical representative, not a Cartesian matrix.
canonicalRotationVectors = zeros(3, numRetained);
orthogonalityErrors = zeros(numRetained, 1);
determinantErrors = zeros(numRetained, 1);
for rotationIndex = 1:numRetained
    Q = canonicalRotations(:,:,rotationIndex);
    logarithm = logm(Q);
    assert(max(abs(imag(logarithm(:)))) < 1e-10, ...
        'The principal rotation logarithm must be real for this representative.');
    logarithm = real(logarithm);
    canonicalRotationVectors(:,rotationIndex) = [ ...
        logarithm(3,2); logarithm(1,3); logarithm(2,1)];
    orthogonalityErrors(rotationIndex) = norm(Q'*Q - eye(n), 'fro');
    determinantErrors(rotationIndex) = abs(det(Q) - 1);
end
canonicalRotationAngles = vecnorm(canonicalRotationVectors, 2, 1).';

assert(isequal(size(canonicalRotations), [n, n, numRetained]));
assert(isequal(size(selectedConditionNumbers), [numRetained, 1]));
assert(isequal(size(selectedIndices), [numRetained, 1]));
assert(all(isfinite(canonicalRotationVectors(:))));
assert(all(diff(selectedConditionNumbers) >= 0));
assert(all(canonicalRotationAngles <= pi + 1e-10));
assert(max(orthogonalityErrors) < 1e-12*max(1,n));
assert(max(determinantErrors) < 1e-12*max(1,n));

fprintf('Selected raw-estimate range: [%.12f, %.12f]\n', ...
    selectedConditionNumbers(1), selectedConditionNumbers(end));
fprintf('Maximum canonical rotation angle: %.12f radians\n', ...
    max(canonicalRotationAngles));

plotLimit = max(0.05, 1.1*max(abs(canonicalRotationVectors(:))));
figure('Name', 'Canonical SO(3) condition-number selection', 'Color', 'w');
scatter3(canonicalRotationVectors(1,:), canonicalRotationVectors(2,:), ...
    canonicalRotationVectors(3,:), 42, selectedConditionNumbers, 'filled', ...
    'DisplayName', '100 canonicalized selected rotations');
hold on;
plot3(0, 0, 0, 'kp', 'MarkerFaceColor', 'k', 'MarkerSize', 10, ...
    'DisplayName', 'Identity');
hold off;
axis equal;
xlim([-plotLimit, plotLimit]);
ylim([-plotLimit, plotLimit]);
zlim([-plotLimit, plotLimit]);
grid on;
view(3);
xlabel('\omega_1 (radians)');
ylabel('\omega_2 (radians)');
zlabel('\omega_3 (radians)');
title('Canonical SO(3) representatives in rotation-vector coordinates');
legend('Location', 'best');
colorbarHandle = colorbar;
ylabel(colorbarHandle, 'Raw pre-canonical \kappa_4 estimate');
