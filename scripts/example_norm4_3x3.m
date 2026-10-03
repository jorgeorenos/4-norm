% Reproducible induced 4-norm estimate and independent spherical comparison in R^3.
scriptDirectory = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDirectory);
addpath(fullfile(projectRoot, 'src', '4-norm'));
addpath(fullfile(projectRoot, 'src', 'helpers'));
rng(42, 'twister');

n = 3;
Q = haar_so(n);
normHandle = @norm_4;
options = struct('NumRandomStarts', 100, 'MaxIterations', 1000);
qNorm4 = compute_induced_norm(Q, normHandle, options);
disp('Matrix Q in SO(3):');
disp(Q);
fprintf('Estimated induced 4-norm: %.12f\n', qNorm4);
% Independent comparison in dimension 3: spherical grid and local refinement.
azimuthCount = 180;
polarCount = 90;
azimuth = (0:azimuthCount) * (2*pi/azimuthCount);
polar = (0:polarCount) * (pi/polarCount);
[azimuthGrid, polarGrid] = meshgrid(azimuth, polar);

sinPolar = sin(polarGrid);
Xgrid = [reshape(cos(azimuthGrid).*sinPolar, 1, []); ...
    reshape(sin(azimuthGrid).*sinPolar, 1, []); ...
    reshape(cos(polarGrid), 1, [])];
Xgrid = Xgrid ./ sum(abs(Xgrid).^4, 1).^(1/4);
amplitudes = reshape(sum(abs(Q*Xgrid).^4, 1).^(1/4), size(azimuthGrid));
[angularReference, bestGridIndex] = max(amplitudes(:));

% Refine the strongest grid directions without using the power iteration.
[~, sortedIndices] = sort(amplitudes(:), 'descend');
numRefinementStarts = min(12, numel(sortedIndices));
refinementOptions = optimset('Display', 'off', 'MaxFunEvals', 1000, ...
    'MaxIter', 1000, 'TolX', 1e-12, 'TolFun', 1e-12);
bestAngles = [azimuthGrid(bestGridIndex); polarGrid(bestGridIndex)];
for j = 1:numRefinementStarts
    index = sortedIndices(j);
    initialAngles = [azimuthGrid(index); polarGrid(index)];
    [refinedAngles, negativeAmplitude] = fminsearch( ...
        @(angles) -spherical_amplitude(angles, Q, normHandle), ...
        initialAngles, refinementOptions);
    if -negativeAmplitude > angularReference
        angularReference = -negativeAmplitude;
        bestAngles = refinedAngles;
    end
end
referenceDirection = [cos(bestAngles(1))*sin(bestAngles(2)); ...
    sin(bestAngles(1))*sin(bestAngles(2)); cos(bestAngles(2))];
xAngularBest = referenceDirection / normHandle(referenceDirection);
fprintf('Approximate spherical reference: %.12f\n', angularReference);
fprintf('Power - spherical difference: %.3e\n', qNorm4-angularReference);

% Visualize the sampled 4-norm sphere and the independently refined maximum.
X1 = reshape(Xgrid(1, :), size(azimuthGrid));
X2 = reshape(Xgrid(2, :), size(azimuthGrid));
X3 = reshape(Xgrid(3, :), size(azimuthGrid));

figure;
surf(X1, X2, X3, amplitudes, 'EdgeColor', 'none', 'FaceColor', 'interp');
hold on;
plot3(xAngularBest(1), xAngularBest(2), xAngularBest(3), 'p', ...
    'Color', [0.85 0.33 0.10], 'MarkerFaceColor', [0.85 0.33 0.10], ...
    'MarkerSize', 12);
axis equal;
grid on;
xlabel('x_1');
ylabel('x_2');
zlabel('x_3');
title(sprintf('4-norm unit sphere | induced norm estimate %.6f', qNorm4));
legend('Spherical sample', 'Refined spherical maximum', 'Location', 'best');
colorbar;
view(3);
hold off;
