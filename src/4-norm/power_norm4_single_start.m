function value = power_norm4_single_start(Q, normHandle, x0, options)
%POWER_NORM4_SINGLE_START Generalized power iteration from one initial direction.
%   The handle must mathematically represent the vector 4-norm; the cubic
%   update is not valid for other vector norms. Returns only the best
%   feasible amplitude (-Inf if no evaluated vector is feasible).
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
    if abs(checked_norm(normHandle, x, 'x', true) - 1) > options.FeasibilityTolerance
        value = -Inf;
    end

    for k = 1:options.MaxIterations
        u = sign(v).*abs(v).^(1/3);
        if any(~isfinite(u)) || ~any(u)
            break;
        end
        uScale = checked_norm(normHandle, u, 'updated u', true);
        xNew = u / uScale;
        if any(~isfinite(xNew))
            break;
        end
        yNew = Q*xNew;
        if any(~isfinite(yNew))
            break;
        end
        valueNew = checked_norm(normHandle, yNew, 'Q*xNew', true);
        feasibilityError = abs(checked_norm(normHandle, xNew, 'xNew', true) - 1);
        if feasibilityError <= options.FeasibilityTolerance && valueNew > value
            value = valueNew;
        end

        if valueNew < oldValue - 128*eps(max([1, abs(valueNew), abs(oldValue)]))
            break;
        end
        vNew = Q'*(yNew.^3);
        lambdaNew = valueNew^4;
        xCube = xNew.^3;
        residual = vNew - lambdaNew*xCube;
        if ~isfinite(lambdaNew) || any(~isfinite(vNew)) || ...
                any(~isfinite(residual))
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
            break;
        end
        v = vNew;
        oldValue = valueNew;
    end

end

function value = checked_norm(normHandle, vector, label, requirePositive)
    value = normHandle(vector);
    if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ...
            ~isfinite(value) || value < 0 || (requirePositive && value == 0)
        error('power_norm4_single_start:InvalidHandleOutput', ...
            'normHandle(%s) must return a finite real scalar, positive for nonzero vectors.', label);
    end
end
