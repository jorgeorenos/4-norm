function values = apply_norm4_handle(normHandle, X, vectorized)
%APPLY_NORM4_HANDLE Evaluate the supplied handle on columns of each page.
%   Only @norm_4 has a known dimension-aware contract. Other equivalent
%   handles are evaluated one vector at a time without probing their API.
    n = size(X, 1);
    if vectorized
        values = normHandle(X, 1);
    else
        columns = reshape(X, n, []);
        values = zeros(1, size(columns, 2));
        for j = 1:size(columns, 2)
            value = normHandle(columns(:,j));
            if ~isnumeric(value) || ~isscalar(value) || ~isreal(value)
                error('compute_induced_norm:InvalidHandleOutput', ...
                    'normHandle must return a finite real scalar for each vector.');
            end
            values(j) = value;
        end
    end
    expectedSize = [1, size(X,2), size(X,3)];
    if ~isnumeric(values) || ~isreal(values) || ...
            numel(values) ~= prod(expectedSize) || any(~isfinite(values(:))) || ...
            any(values(:) < 0)
        error('compute_induced_norm:InvalidHandleOutput', ...
            'normHandle must return finite nonnegative real norms.');
    end
    values = reshape(values, expectedSize);
    if any(reshape(values == 0 & any(X ~= 0, 1), [], 1))
        error('compute_induced_norm:InvalidHandleOutput', ...
            'normHandle must be positive for every nonzero vector.');
    end
end
