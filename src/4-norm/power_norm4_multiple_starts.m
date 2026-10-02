function [values, xHats, runInfo] = power_norm4_multiple_starts(Q, starts, options, includeDiagnostics)
%POWER_NORM4_MULTIPLE_STARTS Vectorized generalized power iterations.
%   This internal kernel is specialized to the 4-norm.  Its columns are
%   independent starts, but their dense matrix products are evaluated in
%   batches.  It is used only when COMPUTE_INDUCED_NORM receives @norm_4.

    if nargin < 4
        includeDiagnostics = true;
    end

    numStarts = size(starts, 2);
    X = starts ./ norm4_columns(starts);
    Y = Q*X;
    oldValues = norm4_columns(Y);
    startValues = oldValues;
    values = oldValues;
    xHats = X;
    candidateFeasibility = abs(norm4_columns(X) - 1);

    % V is the value needed in the following update.  Keeping it avoids
    % repeating Q'*(Y.^3) at the beginning of every iteration.
    V = Q'*(Y.^3);
    converged = false(1, numStarts);
    iterations = zeros(1, numStarts);
    reasonCode = ones(1, numStarts); % 1=max, 2=converged, 3=decrease, 4=failure
    active = true(1, numStarts);
    if options.StoreHistory
        histories = num2cell(oldValues);
    end

    for k = 1:options.MaxIterations
        activeIndices = find(active);
        if isempty(activeIndices)
            break;
        end

        v = V(:, activeIndices);
        u = sign(v).*abs(v).^(1/3);
        validU = all(isfinite(u), 1) & any(u, 1);
        failed = activeIndices(~validU);
        reasonCode(failed) = 4;
        active(failed) = false;

        indices = activeIndices(validU);
        if isempty(indices)
            continue;
        end
        u = u(:, validU);
        xNew = u ./ norm4_columns(u);
        validX = all(isfinite(xNew), 1);
        failed = indices(~validX);
        reasonCode(failed) = 4;
        active(failed) = false;

        indices = indices(validX);
        if isempty(indices)
            continue;
        end
        xNew = xNew(:, validX);
        yNew = Q*xNew;
        valueNew = norm4_columns(yNew);
        validY = all(isfinite(yNew), 1) & isfinite(valueNew) & valueNew > 0;
        failed = indices(~validY);
        reasonCode(failed) = 4;
        active(failed) = false;

        indices = indices(validY);
        if isempty(indices)
            continue;
        end
        xNew = xNew(:, validY);
        yNew = yNew(:, validY);
        valueNew = valueNew(validY);
        feasibilityError = abs(norm4_columns(xNew) - 1);
        improves = feasibilityError <= options.FeasibilityTolerance & ...
            valueNew > values(indices);
        improvedIndices = indices(improves);
        values(improvedIndices) = valueNew(improves);
        xHats(:, improvedIndices) = xNew(:, improves);
        candidateFeasibility(improvedIndices) = feasibilityError(improves);
        iterations(indices) = k;
        if options.StoreHistory
            for j = 1:numel(indices)
                histories{indices(j)}(end+1) = valueNew(j);
            end
        end

        decreaseScale = max([ones(1, numel(indices)); abs(valueNew); ...
            abs(oldValues(indices))], [], 1);
        decreases = valueNew < oldValues(indices) - 128*eps(decreaseScale);
        decreasedIndices = indices(decreases);
        reasonCode(decreasedIndices) = 3;
        active(decreasedIndices) = false;

        indices = indices(~decreases);
        if isempty(indices)
            continue;
        end
        xNew = xNew(:, ~decreases);
        yNew = yNew(:, ~decreases);
        valueNew = valueNew(~decreases);
        feasibilityError = feasibilityError(~decreases);
        vNew = Q'*(yNew.^3);
        lambdaNew = valueNew.^4;
        xCube = xNew.^3;
        residual = vNew - xCube.*lambdaNew;
        validResidual = isfinite(lambdaNew) & all(isfinite(vNew), 1) & ...
            all(isfinite(residual), 1);
        failed = indices(~validResidual);
        reasonCode(failed) = 4;
        active(failed) = false;

        indices = indices(validResidual);
        if isempty(indices)
            continue;
        end
        xNew = xNew(:, validResidual);
        yNew = yNew(:, validResidual);
        valueNew = valueNew(validResidual);
        feasibilityError = feasibilityError(validResidual);
        vNew = vNew(:, validResidual);
        lambdaNew = lambdaNew(validResidual);
        xCube = xCube(:, validResidual);
        residual = residual(:, validResidual);
        denominator = max(norm4_columns(vNew) + abs(lambdaNew).* ...
            norm4_columns(xCube), realmin);
        relativeResidual = norm4_columns(residual) ./ denominator;
        relativeChange = abs(valueNew - oldValues(indices)) ./ ...
            max([ones(1, numel(indices)); abs(valueNew); ...
            abs(oldValues(indices))], [], 1);
        nowConverged = relativeChange <= options.NormTolerance & ...
            relativeResidual <= options.StationarityTolerance & ...
            feasibilityError <= options.FeasibilityTolerance;
        convergedIndices = indices(nowConverged);
        converged(convergedIndices) = true;
        reasonCode(convergedIndices) = 2;
        active(convergedIndices) = false;

        continuing = ~nowConverged;
        continuingIndices = indices(continuing);
        X(:, continuingIndices) = xNew(:, continuing);
        Y(:, continuingIndices) = yNew(:, continuing);
        V(:, continuingIndices) = vNew(:, continuing);
        oldValues(continuingIndices) = valueNew(continuing);
    end

    runInfo.converged = converged;
    runInfo.iterations = iterations;
    runInfo.feasibilityErrors = candidateFeasibility;
    if includeDiagnostics
        runInfo.startValues = startValues;
        runInfo.terminationReasons = reason_labels(reasonCode);
        if options.StoreHistory
            runInfo.histories = histories;
        end
        yBest = Q*xHats;
        vBest = Q'*(yBest.^3);
        lambdaBest = values.^4;
        xCubeBest = xHats.^3;
        residual = vBest - xCubeBest.*lambdaBest;
        denominator = max(norm4_columns(vBest) + abs(lambdaBest).* ...
            norm4_columns(xCubeBest), realmin);
        runInfo.stationarityResiduals = norm4_columns(residual) ./ denominator;
    else
        runInfo.reasonCodes = reasonCode;
    end
end

function values = norm4_columns(X)
    scales = max(abs(X), [], 1);
    values = zeros(1, size(X, 2));
    nonzero = scales > 0;
    if any(nonzero)
        scaled = abs(X(:, nonzero)) ./ scales(nonzero);
        values(nonzero) = scales(nonzero).*sum(scaled.^4, 1).^(1/4);
    end
end

function labels = reason_labels(codes)
    labels = cell(1, numel(codes));
    for j = 1:numel(codes)
        switch codes(j)
            case 1
                labels{j} = 'maxIterations';
            case 2
                labels{j} = 'converged';
            case 3
                labels{j} = 'numericalDecrease';
            otherwise
                labels{j} = 'numericalFailure';
        end
    end
end
