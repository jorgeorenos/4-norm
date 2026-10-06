function qNorm4 = compute_4_norm(Q, normHandle, options)
%COMPUTE_4_NORM Estimate one induced 4-norm per square matrix.
%   Q is a real double n-by-n matrix or n-by-n-by-M array. normHandle must
%   implement the vector 4-norm and return one value per input column and
%   page. The only output is an M-by-1 vector in page order. Every matrix
%   uses independent random starts; no global maximum is certified.
%   MaxWorkingMemoryMB controls the estimated page block size. Screening is
%   an optional approximate two-stage search and can discard a late winner.
    if nargin < 3
        options = struct();
    end
    if ~isa(normHandle, 'function_handle')
        error('compute_4_norm:InvalidHandle', ...
            'normHandle must be a function handle implementing the vector 4-norm.');
    end
    if ~isa(Q, 'double') || ~isreal(Q) || ndims(Q) > 3 || ...
            isempty(Q) || size(Q,1) ~= size(Q,2) || any(~isfinite(Q(:)))
        error('compute_4_norm:InvalidMatrix', ...
            'Q must be a nonempty finite real double n-by-n or n-by-n-by-M array.');
    end
    n = size(Q, 1);
    numMatrices = size(Q, 3);
    defaults = struct('NumRandomStarts', 50, 'MaxIterations', 1000, ...
        'NormTolerance', 1e-10, 'StationarityTolerance', 1e-8, ...
        'FeasibilityTolerance', 1e-12, 'MaxWorkingMemoryMB', 128, ...
        'ScreenIterations', 0, 'NumFinalists', 3);
    if ~isstruct(options) || ~isscalar(options)
        error('compute_4_norm:InvalidOptions', ...
            'options must be a scalar struct.');
    end
    names = fieldnames(options);
    for nameIndex = 1:numel(names)
        name = names{nameIndex};
        if ~isfield(defaults, name)
            error('compute_4_norm:UnknownOption', ...
                'Unknown options field: %s.', name);
        end
        defaults.(name) = options.(name);
    end
    options = defaults;
    if ~valid_integer(options.NumRandomStarts, 1) || ...
            ~valid_integer(options.MaxIterations, 1) || ...
            ~valid_integer(options.ScreenIterations, 0) || ...
            ~valid_integer(options.NumFinalists, 1) || ...
            (options.ScreenIterations > 0 && ...
            (options.NumFinalists > options.NumRandomStarts || ...
            options.ScreenIterations > options.MaxIterations))
        error('compute_4_norm:InvalidCount', ...
            'Invalid start count, iteration count, or screening settings.');
    end
    toleranceNames = {'NormTolerance', 'StationarityTolerance', 'FeasibilityTolerance'};
    for toleranceIndex = 1:numel(toleranceNames)
        tolerance = options.(toleranceNames{toleranceIndex});
        if ~isa(tolerance, 'double') || ~isscalar(tolerance) || ...
                ~isreal(tolerance) || ~isfinite(tolerance) || tolerance <= 0
            error('compute_4_norm:InvalidTolerance', ...
                '%s must be a finite positive double scalar.', ...
                toleranceNames{toleranceIndex});
        end
    end
    if ~isa(options.MaxWorkingMemoryMB, 'double') || ...
            ~isscalar(options.MaxWorkingMemoryMB) || ...
            ~isreal(options.MaxWorkingMemoryMB) || ...
            ~isfinite(options.MaxWorkingMemoryMB) || options.MaxWorkingMemoryMB <= 0
        error('compute_4_norm:InvalidMemoryBudget', ...
            'MaxWorkingMemoryMB must be a finite positive double scalar.');
    end

    bytesPerPage = 8*(n*n + 32*n*options.NumRandomStarts + ...
        16*options.NumRandomStarts);
    blockSize = min(numMatrices, max(1, floor( ...
        options.MaxWorkingMemoryMB*2^20 / bytesPerPage)));
    qNorm4 = zeros(numMatrices, 1);
    if options.ScreenIterations > 0
        screeningOptions = options;
        screeningOptions.MaxIterations = options.ScreenIterations;
        screenedAll = -Inf(options.NumRandomStarts, numMatrices);
        finalistsAll = zeros(n, options.NumFinalists, numMatrices);
        for first = 1:blockSize:numMatrices
            indices = first:min(first + blockSize - 1, numMatrices);
            starts = random_starts(n, options.NumRandomStarts, numel(indices));
            screened = power_norm4_multiple_starts( ...
                Q(:,:,indices), normHandle, starts, screeningOptions);
            screened(~isfinite(screened)) = -Inf;
            screenedAll(:,indices) = screened;
            [~, ranking] = sort(screened, 1, 'descend');
            flatStarts = reshape(starts, n, []);
            linearIds = ranking(1:options.NumFinalists,:) + ...
                (0:numel(indices)-1)*options.NumRandomStarts;
            finalistsAll(:,:,indices) = reshape(flatStarts(:,linearIds), ...
                n, options.NumFinalists, numel(indices));
        end
        finalistBytesPerPage = 8*(n*n + 32*n*options.NumFinalists + ...
            16*options.NumFinalists);
        finalistBlockSize = min(numMatrices, max(1, floor( ...
            options.MaxWorkingMemoryMB*2^20 / finalistBytesPerPage)));
        refinedAll = -Inf(options.NumFinalists, numMatrices);
        for first = 1:finalistBlockSize:numMatrices
            indices = first:min(first + finalistBlockSize - 1, numMatrices);
            refined = power_norm4_multiple_starts( ...
                Q(:,:,indices), normHandle, finalistsAll(:,:,indices), options);
            refined(~isfinite(refined)) = -Inf;
            refinedAll(:,indices) = refined;
        end
        estimates = max(max(screenedAll, [], 1), max(refinedAll, [], 1));
        failed = find(~isfinite(estimates), 1);
        if ~isempty(failed)
            error('compute_4_norm:NoCandidate', ...
                'No start produced a finite feasible candidate for Q(:,:,%d).', failed);
        end
        qNorm4 = estimates(:);
        return;
    end

    for first = 1:blockSize:numMatrices
        indices = first:min(first + blockSize - 1, numMatrices);
        starts = random_starts(n, options.NumRandomStarts, numel(indices));
        candidates = power_norm4_multiple_starts( ...
            Q(:,:,indices), normHandle, starts, options);
        candidates(~isfinite(candidates)) = -Inf;
        estimates = max(candidates, [], 1);
        failed = find(~isfinite(estimates), 1);
        if ~isempty(failed)
            error('compute_4_norm:NoCandidate', ...
                'No start produced a finite feasible candidate for Q(:,:,%d).', ...
                indices(failed));
        end
        qNorm4(indices) = estimates(:);
    end
end

function starts = random_starts(n, numStarts, numMatrices)
    starts = randn(n, numStarts, numMatrices);
    columns = reshape(starts, n, []);
    zeroColumns = ~any(columns, 1);
    while any(zeroColumns)
        columns(:,zeroColumns) = randn(n, sum(zeroColumns));
        zeroColumns = ~any(columns, 1);
    end
    starts = reshape(columns, n, numStarts, numMatrices);
end

function valid = valid_integer(value, minimum)
    valid = isa(value, 'double') && isscalar(value) && isreal(value) && ...
        isfinite(value) && value >= minimum && value == floor(value);
end
