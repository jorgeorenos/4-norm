function values = norm4_columns(X)
%NORM4_COLUMNS Evaluate scaled 4-norms along the first dimension.
%   X is an n-by-S array or n-by-S-by-B array. This function can be passed
%   directly as the normHandle argument of the handle-based 4-norm routines.
%   The direct path is used for ordinary magnitudes; a scaled fallback avoids
%   unnecessary overflow and underflow.
    squares = X.*X;
    values = sum(squares.*squares, 1).^(1/4);
    if all(isfinite(values(:))) && all(values(:) > 0)
        return;
    end
    nonzero = any(X ~= 0, 1);
    if all(isfinite(values(:))) && all((values > 0 | ~nonzero), 'all')
        return;
    end
    magnitudes = abs(X);
    scales = max(magnitudes, [], 1);
    divisors = scales;
    divisors(divisors == 0) = 1;
    scaled = magnitudes ./ divisors;
    scaledSquares = scaled.*scaled;
    values = scales .* sum(scaledSquares.*scaledSquares, 1).^(1/4);
end
