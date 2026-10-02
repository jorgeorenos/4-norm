function X = generate_unit_l4_vectors(n, N, method)
%GENERATE_UNIT_L4_VECTORS Generate columns on the 4-norm unit sphere.
%   X = GENERATE_UNIT_L4_VECTORS(N, M) returns M vectors of dimension N
%   with unit 4-norm. The random method normalizes Gaussian directions
%   and produces valid points, but is not guaranteed to be uniform with
%   respect to surface area or cone measure.
%
%   X = GENERATE_UNIT_L4_VECTORS(N, M, METHOD) accepts METHOD as
%   'random' (default) or 'angular'. The 'angular' method is only
%   available in dimension 2 and orders equally spaced directions.

    if nargin < 3
        method = 'random';
    end

    validate_positive_integer(n, 'n');
    validate_positive_integer(N, 'N');
    if ~ischar(method) || ~ismember(method, {'random', 'angular'})
        error('generate_unit_l4_vectors:InvalidMethod', ...
            "method must be 'random' or 'angular'.");
    end

    switch method
        case 'random'
            Z = randn(n, N);
            scales = sum(abs(Z).^4, 1).^(1/4);
            while any(scales == 0)
                zeroColumns = (scales == 0);
                Z(:, zeroColumns) = randn(n, sum(zeroColumns));
                scales = sum(abs(Z).^4, 1).^(1/4);
            end
            X = Z ./ scales;

        case 'angular'
            if n ~= 2
                error('generate_unit_l4_vectors:AngularDimension', ...
                    "The 'angular' method is only available when n is 2.");
            end
            if N < 4
                error('generate_unit_l4_vectors:AngularCount', ...
                    "The 'angular' method requires N to be at least 4.");
            end
            theta = (0:N-1) * (2*pi/N);
            Z = [cos(theta); sin(theta)];
            X = Z ./ sum(abs(Z).^4, 1).^(1/4);
    end
end

function validate_positive_integer(value, name)
    if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ...
            ~isfinite(value) || value <= 0 || value ~= floor(value)
        error('generate_unit_l4_vectors:InvalidSize', ...
            '%s must be a finite positive integer scalar.', name);
    end
end
