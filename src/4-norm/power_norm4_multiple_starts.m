function values = power_norm4_multiple_starts(Q, starts, options)
%POWER_NORM4_MULTIPLE_STARTS Vectorized generalized power iterations.
%   This internal kernel is specialized to the 4-norm.  Its columns are
%   independent starts, but their dense matrix products are evaluated in
%   batches.  It is used only when COMPUTE_INDUCED_NORM receives @norm_4.
%   Returns one best feasible amplitude per start (-Inf if none is feasible).

    numStarts = size(starts, 2);
    X = starts ./ norm4_columns(starts);
    Y = Q*X;
    oldValues = norm4_columns(Y);
    values = oldValues;
    values(abs(norm4_columns(X) - 1) > options.FeasibilityTolerance) = -Inf;

    % V is the value needed in the following update.  Keeping it avoids
    % repeating Q'*(Y.^3) at the beginning of every iteration.
    V = Q'*(Y.^3);
    active = true(1, numStarts);

    for k = 1:options.MaxIterations
        activeIndices = find(active);
        if isempty(activeIndices)
            break;
        end

        v = V(:, activeIndices);
        u = sign(v).*abs(v).^(1/3);
        validU = all(isfinite(u), 1) & any(u, 1);
        failed = activeIndices(~validU);
        active(failed) = false;

        indices = activeIndices(validU);
        if isempty(indices)
            continue;
        end
        u = u(:, validU);
        xNew = u ./ norm4_columns(u);
        validX = all(isfinite(xNew), 1);
        failed = indices(~validX);
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

        decreaseScale = max([ones(1, numel(indices)); abs(valueNew); ...
            abs(oldValues(indices))], [], 1);
        decreases = valueNew < oldValues(indices) - 128*eps(decreaseScale);
        decreasedIndices = indices(decreases);
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
        active(failed) = false;

        indices = indices(validResidual);
        if isempty(indices)
            continue;
        end
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
        active(convergedIndices) = false;

        continuing = ~nowConverged;
        continuingIndices = indices(continuing);
        V(:, continuingIndices) = vNew(:, continuing);
        oldValues(continuingIndices) = valueNew(continuing);
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
