function qNorm4 = compute_induced_norm(Q, normHandle, options)
%COMPUTE_INDUCED_NORM Estimate one induced 4-norm per matrix in SO(n).
%   Q is a real double n-by-n matrix or n-by-n-by-M array. The only
%   output is an M-by-1 vector in the same order as Q(:,:,k).
%   normHandle must mathematically implement the vector 4-norm.
%   Each matrix uses independent random starts; no global maximum is certified.
    if nargin < 3
        options = struct();
    end
    if ~isa(normHandle, 'function_handle')
        error('compute_induced_norm:InvalidHandle', ...
            'normHandle must be a function handle for the 4-norm.');
    end
    if ~isa(Q, 'double') || ~isreal(Q) || ndims(Q) > 3 || ...
            isempty(Q) || size(Q,1) ~= size(Q,2) || any(~isfinite(Q(:)))
        error('compute_induced_norm:InvalidMatrix', ...
            'Q must be a nonempty finite real double n-by-n or n-by-n-by-M array.');
    end
    n = size(Q, 1);
    numMatrices = size(Q, 3);
    matrixTolerance = 1e-12*max(1,n);
    % Validate the whole batch before consuming the random state.
    for k = 1:numMatrices
        matrix = Q(:,:,k);
        if norm(matrix'*matrix-eye(n), 'fro') > matrixTolerance || ...
                abs(det(matrix)-1) > matrixTolerance
            error('compute_induced_norm:NotSO', ...
                'Q(:,:,%d) must belong to SO(n) within tolerance %.3e.', ...
                k, matrixTolerance);
        end
    end
    defaults = struct('NumRandomStarts', 50, 'MaxIterations', 1000, ...
        'NormTolerance', 1e-10, 'StationarityTolerance', 1e-8, ...
        'FeasibilityTolerance', 1e-12);
    if ~isstruct(options) || ~isscalar(options)
        error('compute_induced_norm:InvalidOptions', ...
            'options must be a scalar struct.');
    end
    names = fieldnames(options);
    for j = 1:numel(names)
        name = names{j};
        if ~isfield(defaults, name)
            error('compute_induced_norm:UnknownOption', ...
                'Unknown options field: %s.', name);
        end
        defaults.(name) = options.(name);
    end
    options = defaults;
    if ~valid_integer(options.NumRandomStarts, 1) || ...
            ~valid_integer(options.MaxIterations, 1)
        error('compute_induced_norm:InvalidCount', ...
            'NumRandomStarts and MaxIterations must be positive integers.');
    end
    tolerances = {'NormTolerance', 'StationarityTolerance', 'FeasibilityTolerance'};
    for j = 1:numel(tolerances)
        tolerance = options.(tolerances{j});
        if ~isa(tolerance, 'double') || ~isscalar(tolerance) || ...
                ~isreal(tolerance) || ~isfinite(tolerance) || tolerance <= 0
            error('compute_induced_norm:InvalidTolerance', ...
                '%s must be a finite positive double scalar.', tolerances{j});
        end
    end

    qNorm4 = zeros(numMatrices, 1);
    useVectorizedKernel = strcmp(func2str(normHandle), 'norm_4');
    for k = 1:numMatrices
        matrix = Q(:,:,k);
        starts = generate_unit_l4_vectors(n, options.NumRandomStarts, 'random');
        for j = 1:options.NumRandomStarts
            starts(:,j) = starts(:,j) / checked_norm(normHandle, starts(:,j), 'random start');
        end
        if useVectorizedKernel
            candidates = power_norm4_multiple_starts(matrix, starts, options);
        else
            candidates = zeros(1, options.NumRandomStarts);
            for j = 1:options.NumRandomStarts
                candidates(j) = power_norm4_single_start( ...
                    matrix, normHandle, starts(:,j), options);
            end
        end
        feasibleValues = candidates(isfinite(candidates));
        if isempty(feasibleValues)
            error('compute_induced_norm:NoCandidate', ...
                'No start produced a finite feasible candidate for Q(:,:,%d).', k);
        end
        qNorm4(k) = max(feasibleValues);
        if qNorm4(k) < 1-matrixTolerance || qNorm4(k) > n^(1/4)+matrixTolerance
            warning('compute_induced_norm:BoundAnomaly', ...
                'The estimate for Q(:,:,%d) falls outside the theoretical bounds.', k);
        end
    end
end

function valid = valid_integer(value, minimum)
    valid = isa(value, 'double') && isscalar(value) && isreal(value) && ...
        isfinite(value) && value >= minimum && value == floor(value);
end

function value = checked_norm(normHandle, vector, label)
    value = normHandle(vector);
    if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ...
            ~isfinite(value) || value <= 0
        error('compute_induced_norm:InvalidHandleOutput', ...
            'normHandle(%s) must return a finite real scalar, positive for nonzero vectors.', label);
    end
end
