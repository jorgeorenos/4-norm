function [value, xHat, runInfo] = power_norm4_single_start(Q, normHandle, x0, options)
%POWER_NORM4_SINGLE_START Generalized power iteration from one initial direction.
%   The handle must mathematically represent the vector 4-norm; the cubic
%   update is not valid for other vector norms.
    if ~isa(normHandle, 'function_handle')
        error('power_norm4_single_start:InvalidHandle', ...
            'normHandle must be a function handle for the 4-norm.');
    end
    n = size(Q, 1);
    if ~isa(x0, 'double') || ~isvector(x0) || numel(x0) ~= n || ...
            ~isreal(x0) || any(~isfinite(x0(:)))
        error('power_norm4_single_start:InvalidStart', ...
            'x0 must be a finite real double vector of dimension n.');
    end
    x = x0(:);
    xScale = checked_norm(normHandle, x, 'x0', true);
    x = x / xScale;
    y = Q*x;
    oldValue = checked_norm(normHandle, y, 'initial Q*x', true);
    % Keep the gradient-like quantity for the next iteration.  The old
    % implementation recomputed both y and v immediately after having
    % obtained them as yNew and vNew in the previous iteration.
    v = Q'*(y.^3);
    value = oldValue;
    xHat = x;
    runInfo.startValue = oldValue;
    runInfo.converged = false;
    runInfo.terminationReason = 'maxIterations';
    runInfo.iterations = 0;
    if options.StoreHistory
        runInfo.history = oldValue;
    end

    for k = 1:options.MaxIterations
        u = sign(v).*abs(v).^(1/3);
        if any(~isfinite(u)) || ~any(u)
            runInfo.terminationReason = 'numericalFailure';
            break;
        end
        uScale = checked_norm(normHandle, u, 'updated u', true);
        xNew = u / uScale;
        if any(~isfinite(xNew))
            runInfo.terminationReason = 'numericalFailure';
            break;
        end
        yNew = Q*xNew;
        if any(~isfinite(yNew))
            runInfo.terminationReason = 'numericalFailure';
            break;
        end
        valueNew = checked_norm(normHandle, yNew, 'Q*xNew', true);
        feasibilityError = abs(checked_norm(normHandle, xNew, 'xNew', true) - 1);
        if feasibilityError <= options.FeasibilityTolerance && valueNew > value
            value = valueNew;
            xHat = xNew;
        end
        runInfo.iterations = k;
        if options.StoreHistory
            runInfo.history(end+1) = valueNew;
        end
        if valueNew < oldValue - 128*eps(max([1, abs(valueNew), abs(oldValue)]))
            runInfo.terminationReason = 'numericalDecrease';
            break;
        end
        vNew = Q'*(yNew.^3);
        lambdaNew = valueNew^4;
        xCube = xNew.^3;
        residual = vNew - lambdaNew*xCube;
        if ~isfinite(lambdaNew) || any(~isfinite(vNew)) || ...
                any(~isfinite(residual))
            runInfo.terminationReason = 'numericalFailure';
            break;
        end
        denominator = max(checked_norm(normHandle, vNew, 'vNew', false) + ...
            abs(lambdaNew)*checked_norm(normHandle, xCube, 'xNew.^3', true), realmin);
        relativeResidual = checked_norm(normHandle, residual, 'residual', false) / denominator;
        relativeChange = abs(valueNew - oldValue) / ...
            max([1, abs(valueNew), abs(oldValue)]);
        if relativeChange <= options.NormTolerance && ...
                relativeResidual <= options.StationarityTolerance && ...
                feasibilityError <= options.FeasibilityTolerance
            runInfo.converged = true;
            runInfo.terminationReason = 'converged';
            break;
        end
        x = xNew;
        y = yNew;
        v = vNew;
        oldValue = valueNew;
    end

    % All diagnostics correspond to the candidate actually returned.
    yBest = Q*xHat;
    vBest = Q'*(yBest.^3);
    lambdaBest = value^4;
    xCubeBest = xHat.^3;
    if all(isfinite(vBest)) && isfinite(lambdaBest)
        bestResidual = vBest - lambdaBest*xCubeBest;
        denominator = max(checked_norm(normHandle, vBest, 'vBest', false) + ...
            abs(lambdaBest)*checked_norm(normHandle, xCubeBest, 'xHat.^3', true), realmin);
        runInfo.stationarityResidual = ...
            checked_norm(normHandle, bestResidual, 'final residual', false) / denominator;
    else
        runInfo.stationarityResidual = Inf;
    end
    runInfo.feasibilityError = abs(checked_norm(normHandle, xHat, 'xHat', true) - 1);
end

function value = checked_norm(normHandle, vector, label, requirePositive)
    value = normHandle(vector);
    if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ...
            ~isfinite(value) || value < 0 || (requirePositive && value == 0)
        error('power_norm4_single_start:InvalidHandleOutput', ...
            'normHandle(%s) must return a finite real scalar, positive for nonzero vectors.', label);
    end
end
