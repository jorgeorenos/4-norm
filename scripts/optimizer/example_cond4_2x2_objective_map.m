% Explore low condition-number representatives in the canonical SO(2) domain.
% A dense angular reference is compared with low estimated condition-number
% rotations selected from a Haar sweep.
scriptDirectory = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptDirectory));
addpath(fullfile(projectRoot, 'src', '4-norm_handle'), '-begin');
addpath(fullfile(projectRoot, 'src', 'helpers'), '-begin');

numTheta = 181;
numPhi = 1440;
numSamples = 10000;
numRetained = 100;
batchSize = 100;
numRandomStarts = 75;
options = struct('NumRandomStarts', numRandomStarts, 'MaxIterations', 1000, ...
    'MaxWorkingMemoryMB', 128, 'ScreenIterations', 2, 'NumFinalists', 4);
normHandle = @norm4_columns;

canonicalTheta = linspace(-pi/4, pi/4, numTheta);
phi = (0:numPhi-1)*(2*pi/numPhi);
directions = [cos(phi); sin(phi)];
X = directions ./ norm4_columns(directions);

referenceConditionNumbers = zeros(numTheta, 1);
for thetaIndex = 1:numTheta
    angle = canonicalTheta(thetaIndex);
    Q = [cos(angle), -sin(angle); sin(angle), cos(angle)];
    referenceConditionNumbers(thetaIndex) = ...
        max(norm4_columns(Q*X)) * max(norm4_columns(Q'*X));
end
zeroIndex = find(abs(canonicalTheta) < 16*eps, 1);

rng(42, 'twister');
conditionNumberFn = @(rotations) evaluate_condition_number_batch( ...
    rotations, normHandle, options);
[canonicalRotations, selectedConditionNumbers, selectedIndices] = ...
    initialize_condition_number(2, numSamples, batchSize, ...
    conditionNumberFn, numRetained);
canonicalAngles = reshape(atan2(canonicalRotations(2,1,:), ...
    canonicalRotations(1,1,:)), [], 1);

assert(abs(referenceConditionNumbers(zeroIndex) - 1) < 1e-12);
assert(abs(referenceConditionNumbers(1) - sqrt(2)) < 1e-12);
assert(isequal(size(canonicalRotations), [2, 2, numRetained]));
assert(isequal(size(selectedConditionNumbers), [numRetained, 1]));
assert(isequal(size(selectedIndices), [numRetained, 1]));
assert(all(diff(selectedConditionNumbers) >= 0));
assert(all(abs(canonicalAngles) <= pi/4 + 1e-12));

[minimumCondition, minimumIndex] = min(referenceConditionNumbers);
[maximumCondition, maximumIndex] = max(referenceConditionNumbers);
fprintf('Minimum angular-grid kappa_4: %.12f at theta = %.12f radians\n', ...
    minimumCondition, canonicalTheta(minimumIndex));
fprintf('Maximum angular-grid kappa_4: %.12f at theta = %.12f radians\n', ...
    maximumCondition, canonicalTheta(maximumIndex));
fprintf('Selected raw-estimate range: [%.12f, %.12f]\n', ...
    selectedConditionNumbers(1), selectedConditionNumbers(end));

figure('Name', 'Canonical SO(2) condition-number selection', 'Color', 'w');
plot(canonicalTheta, referenceConditionNumbers, 'LineWidth', 1.5, ...
    'Color', [0.85, 0.40, 0.10], ...
    'DisplayName', 'Angular-grid reference');
hold on;
scatter(canonicalAngles, selectedConditionNumbers, 24, ...
    [0.15, 0.45, 0.75], 'filled', ...
    'DisplayName', '100 canonicalized selected rotations');
plot(canonicalTheta(minimumIndex), minimumCondition, 'ko', ...
    'MarkerFaceColor', 'k', 'DisplayName', 'Reference minimum');
hold off;
grid on;
xlim([-pi/4, pi/4]);
ylim([0.99, sqrt(2) + 0.01]);
set(gca, 'XTick', [-pi/4, 0, pi/4], ...
    'XTickLabel', {'-\pi/4', '0', '\pi/4'});
xlabel('Canonical rotation angle \theta');
ylabel('\kappa_4(Q(\theta))');
title('Selected rotations in the canonical SO(2) domain');
legend('Location', 'northwest');
