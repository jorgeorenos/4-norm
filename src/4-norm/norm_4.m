function value = norm_4(x)
%NORM_4 Vector 4-norm with scaling to avoid unnecessary extremes.
%   Accepts a finite real double row or column vector, including zero.
    if ~isa(x, 'double') || ~isvector(x) || isempty(x) || ...
            ~isreal(x) || any(~isfinite(x(:)))
        error('norm_4:InvalidVector', ...
            'x must be a nonempty, finite, real double vector.');
    end
    s = max(abs(x(:)));
    if s == 0
        value = 0;
    else
        value = s * sum((abs(x(:)) / s).^4)^(1/4);
    end
end
