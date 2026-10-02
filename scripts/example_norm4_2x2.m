% Reproducible induced 4-norm estimate and angular comparison in R^2.
scriptDirectory = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDirectory);
addpath(fullfile(projectRoot, 'src', '4-norm'));
addpath(fullfile(projectRoot, 'src', 'helpers'));
rng(42, 'twister');

n = 2;
Q = haar_so(n);
normHandle = @norm_4;
options = struct('NumRandomStarts', 50, 'MaxIterations', 1000, ...
    'StoreHistory', true);
[qNorm4, xBest, info] = compute_induced_norm(Q, normHandle, options);
disp('Matrix Q in SO(2):');
disp(Q);
fprintf('Estimated induced 4-norm: %.12f\n', qNorm4);
fprintf('4-norm of the candidate vector: %.12f\n', normHandle(xBest));
fprintf('Amplitude of Q*xBest: %.12f\n', normHandle(Q*xBest));
disp('Candidate vector:');
disp(xBest);
disp(info);

% Independent comparison in dimension 2: grid and periodic refinement.
N = 10000;
theta = (0:N-1) * (2*pi/N);
X = generate_unit_l4_vectors(2, N, 'angular');
amplitudes = zeros(1,N);
for k = 1:N
    x = X(:,k) / normHandle(X(:,k));
    amplitudes(k) = normHandle(Q*x);
end
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

bestAngle = mod(atan2(xBest(2), xBest(1)), 2*pi);
figure;
plot([theta, 2*pi], [amplitudes, amplitudes(1)], ...
    'Color', [0.15 0.35 0.70], 'LineWidth', 1.5);
hold on;
plot(bestAngle, qNorm4, 'o', 'Color', [0.85 0.33 0.10], ...
    'MarkerFaceColor', [0.85 0.33 0.10], 'MarkerSize', 7);
xlim([0 2*pi]);
grid on;
xlabel('Angle (radians)');
ylabel('||Qx(θ)||_4');
title('Angular amplitude and best power-method candidate');
legend('Angular sample', 'Power-method candidate', 'Location', 'best');
hold off;
