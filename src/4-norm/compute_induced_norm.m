function [qNorm4, xBest, info] = compute_induced_norm(Q, normHandle, options)
%COMPUTE_INDUCED_NORM Estimate the induced 4-norm of Q in SO(n).
%   qNorm4 = COMPUTE_INDUCED_NORM(Q, @norm_4) returns an estimate.
%   [qNorm4, xBest, info] = COMPUTE_INDUCED_NORM(Q, normHandle, options)
%   provides a unit candidate and per-start diagnostics. normHandle must
%   mathematically implement the 4-norm; other norms are not optimized.
    if nargin < 3
        options = struct();
    end
    if ~isa(normHandle, 'function_handle')
        error('compute_induced_norm:InvalidHandle', ...
            'normHandle must be a function handle for the 4-norm.');
    end
    if ~isa(Q, 'double') || ~isreal(Q) || ~ismatrix(Q) || ...
            isempty(Q) || size(Q,1) ~= size(Q,2) || any(~isfinite(Q(:)))
        error('compute_induced_norm:InvalidMatrix', ...
            'Q must be a nonempty, finite, real square double matrix.');
    end
    n = size(Q, 1);
    matrixTolerance = 1e-12*max(1,n);
    if norm(Q'*Q-eye(n), 'fro') > matrixTolerance || ...
            abs(det(Q)-1) > matrixTolerance
        error('compute_induced_norm:NotSO', ...
            'Q must belong to SO(n) within tolerance %.3e.', matrixTolerance);
    end
    defaults = struct('NumRandomStarts', 50, 'MaxIterations', 1000, ...
        'NormTolerance', 1e-10, 'StationarityTolerance', 1e-8, ...
        'FeasibilityTolerance', 1e-12, 'StoreHistory', false);
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
    if ~valid_integer(options.NumRandomStarts, 0) || ...
            ~valid_integer(options.MaxIterations, 1)
        error('compute_induced_norm:InvalidCount', ...
            'NumRandomStarts must be a nonnegative integer and MaxIterations a positive integer.');
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
    if ~islogical(options.StoreHistory) || ~isscalar(options.StoreHistory)
        error('compute_induced_norm:InvalidHistory', ...
            'StoreHistory must be a logical scalar.');
    end

    numStarts = n + options.NumRandomStarts;
    starts = zeros(n, numStarts);
    for j = 1:n
        z = Q(j,:).';
        scale = checked_norm(normHandle, z, 'structured start');
        starts(:,j) = z / scale;
    end
    if options.NumRandomStarts > 0
        randomStarts = generate_unit_l4_vectors(n, options.NumRandomStarts, 'random');
        for j = 1:options.NumRandomStarts
            z = randomStarts(:,j);
            starts(:,n+j) = z / checked_norm(normHandle, z, 'random start');
        end
    end

    info.method = 'generalized power method with multiple starts';
    info.isEstimate = true;
    info.numStarts = numStarts;
    info.numRandomStarts = options.NumRandomStarts;
    info.startValues = NaN(1,numStarts);
    info.startConverged = false(1,numStarts);
    info.startTerminationReasons = cell(1,numStarts);
    if options.StoreHistory
        info.histories = cell(1,numStarts);
    end
    info.theoreticalLowerBound = 1;
    info.theoreticalUpperBound = n^(1/4);
    qNorm4 = -Inf;
    xBest = [];
    for j = 1:numStarts
        [candidate, x, runInfo] = power_norm4_single_start( ...
            Q, normHandle, starts(:,j), options);
        info.startValues(j) = runInfo.startValue;
        info.startConverged(j) = runInfo.converged;
        info.startTerminationReasons{j} = runInfo.terminationReason;
        if options.StoreHistory
            info.histories{j} = runInfo.history;
        end
        if isfinite(candidate) && isfinite(runInfo.feasibilityError) && ...
                runInfo.feasibilityError <= options.FeasibilityTolerance && ...
                candidate > qNorm4
            qNorm4 = candidate;
            xBest = x;
            info.bestStartIndex = j;
            info.bestConverged = runInfo.converged;
            info.bestTerminationReason = runInfo.terminationReason;
            info.bestIterations = runInfo.iterations;
            info.bestStationarityResidual = runInfo.stationarityResidual;
            info.bestFeasibilityError = runInfo.feasibilityError;
        end
    end
    if isempty(xBest)
        error('compute_induced_norm:NoCandidate', ...
            'No start produced a finite feasible candidate.');
    end
    actualValue = checked_norm(normHandle, Q*xBest, 'Q*xBest');
    if abs(actualValue-qNorm4) > 128*eps(max(1,actualValue))
        error('compute_induced_norm:InconsistentResult', ...
            'The amplitude of the best vector does not match qNorm4.');
    end
    boundTolerance = 1e-12*max(1,n);
    if qNorm4 < info.theoreticalLowerBound-boundTolerance || ...
            qNorm4 > info.theoreticalUpperBound+boundTolerance
        warning('compute_induced_norm:BoundAnomaly', ...
            'The estimate falls outside the theoretical bounds for SO(n).');
    end
    if ~info.bestConverged
        warning('compute_induced_norm:NotConverged', ...
            'The best candidate stopped due to %s after %d iterations.', ...
            info.bestTerminationReason, info.bestIterations);
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
