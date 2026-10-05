% Reproducible induced 4-norm estimate and angular comparison in R^2.
scriptDirectory = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDirectory);
addpath(fullfile(projectRoot, 'src', '4-norm'));
addpath(fullfile(projectRoot, 'src', 'helpers'));
rng(42, 'twister');

n = 2;
Q = haar_so(n);
normHandle = @norm_4;
options = struct('NumRandomStarts', 50, 'MaxIterations', 1000);
qNorm4 = compute_induced_norm(Q, normHandle, options);
disp('Matrix Q in SO(2):');
disp(Q);
fprintf('Estimated induced 4-norm: %.12f\n', qNorm4);
% Independent comparison in dimension 2: grid and periodic refinement.
N = 10000;
theta = (0:N-1) * (2*pi/N);
directions = [cos(theta); sin(theta)];
X = directions ./ normHandle(directions, 1);
% Evaluate the independent angular grid in one matrix product and
% column-wise norm calls instead of 10000 scalar handle invocations.
amplitudes = normHandle(Q*X, 1);
[angularReference, bestGridIndex] = max(amplitudes);
left = amplitudes([N, 1:N-1]);
right = amplitudes([2:N, 1]);
localIndices = find(amplitudes >= left & amplitudes >= right & ...
    (amplitudes > left | amplitudes > right));
if isempty(localIndices)
    localIndices = bestGridIndex; % Flat grid or complete tie.
end
angleStep = 2*pi/N;
for k = localIndices
    [refinedAngle, negativeAmplitude] = fminbnd( ...
        @(angle) -angular_amplitude(angle, Q, normHandle), ...
        theta(k)-angleStep, theta(k)+angleStep);
    if -negativeAmplitude > angularReference
        angularReference = -negativeAmplitude;
    end
end
fprintf('Approximate angular reference: %.12f\n', angularReference);
fprintf('Power - angular difference: %.3e\n', qNorm4-angularReference);

figure;
plot([theta, 2*pi], [amplitudes, amplitudes(1)], ...
    'Color', [0.15 0.35 0.70], 'LineWidth', 1.5);
hold on;
yline(qNorm4, '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.5);
xlim([0 2*pi]);
grid on;
xlabel('Angle (radians)');
ylabel('||Qx(θ)||_4');
title('Angular amplitude and induced 4-norm estimate');
legend('Angular sample', 'Induced 4-norm estimate', 'Location', 'best');
hold off;
