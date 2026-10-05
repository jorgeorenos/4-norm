function value = norm_4(x, dimension)
%NORM_4 Vector 4-norm with scaling to avoid unnecessary extremes.
%   NORM_4(x) accepts a finite real double row or column vector.
%   NORM_4(X, dimension) evaluates vectors along that dimension, including
%   all columns of each page with NORM_4(X, 1). Zero vectors return zero.
    if ~isa(x, 'double') || isempty(x) || ...
            ~isreal(x) || any(~isfinite(x(:)))
        error('norm_4:InvalidVector', ...
            'x must contain nonempty, finite, real double vectors.');
    end
    if nargin < 2
        if ~isvector(x)
            error('norm_4:InvalidVector', ...
                'Specify a dimension to evaluate a matrix or an array.');
        end
        x = x(:);
        dimension = 1;
    elseif ~isnumeric(dimension) || ~isscalar(dimension) || ...
            ~isreal(dimension) || ~isfinite(dimension) || ...
            dimension < 1 || dimension ~= floor(dimension)
        error('norm_4:InvalidDimension', ...
            'dimension must be a finite positive integer scalar.');
    end
    magnitudes = abs(x);
    scales = max(magnitudes, [], dimension);
    divisors = scales;
    divisors(divisors == 0) = 1;
    value = scales .* sum((magnitudes ./ divisors).^4, dimension).^(1/4);
end
