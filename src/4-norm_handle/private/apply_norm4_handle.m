function values = apply_norm4_handle(normHandle, X)
%APPLY_NORM4_HANDLE Evaluate a vector 4-norm handle along columns and pages.
%   The preferred contract accepts an n-by-S-by-B input and returns a
%   1-by-S-by-B array. Handles that accept only one vector are supported by
%   a column-by-column fallback. The handle must implement the vector 4-norm.
    expectedSize = [1, size(X,2), size(X,3)];
    try
        values = normHandle(X);
    catch
        columns = reshape(X, size(X,1), []);
        values = zeros(1, size(columns,2));
        for columnIndex = 1:size(columns,2)
            value = normHandle(columns(:,columnIndex));
            if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ...
                    ~isfinite(value) || value < 0
                error('compute_4_norm:InvalidHandleOutput', ...
                    ['normHandle must return one finite, real, nonnegative ', ...
                    'scalar per nonzero vector.']);
            end
            values(columnIndex) = double(value);
        end
        values = reshape(values, expectedSize);
    end

    if ~isnumeric(values) || ~isreal(values) || ...
            numel(values) ~= prod(expectedSize) || any(~isfinite(values(:))) || ...
            any(values(:) < 0)
        error('compute_4_norm:InvalidHandleOutput', ...
            ['normHandle must return finite, real, nonnegative norms with ', ...
            'one value per input column and page.']);
    end
    values = reshape(double(values), expectedSize);
    nonzero = any(X ~= 0, 1);
    if any(values == 0 & nonzero, 'all')
        error('compute_4_norm:InvalidHandleOutput', ...
            'normHandle must be positive for every nonzero vector.');
    end
end
