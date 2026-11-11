function estimates = compute_norm4_fixed_starts(A, normHandle, starts, options)
%COMPUTE_NORM4_FIXED_STARTS Estimate induced 4-norms from common starts.
    if ~isa(A, 'double') || ~isreal(A) || ndims(A) > 3 || isempty(A) || ...
            size(A,1) ~= size(A,2) || any(~isfinite(A(:)))
        error('optimizer:compute_norm4_fixed_starts:InvalidMatrix', ...
            'A must be a nonempty finite real double square matrix or page array.');
    end
    if ~isa(normHandle, 'function_handle')
        error('optimizer:compute_norm4_fixed_starts:InvalidHandle', ...
            'normHandle must be a function handle.');
    end
    n = size(A, 1);
    if ~isa(starts, 'double') || ~isreal(starts) || ...
            ~ismatrix(starts) || size(starts,1) ~= n || isempty(starts) || ...
            any(~isfinite(starts(:))) || any(~any(starts ~= 0, 1))
        error('optimizer:compute_norm4_fixed_starts:InvalidStarts', ...
            'starts must be a finite real double n-by-S matrix with nonzero columns.');
    end
    options = complete_options(options, size(starts,2));

    numMatrices = size(A, 3);
    numStarts = size(starts, 2);
    bytesPerPage = 8*(n*n + 32*n*numStarts + 16*numStarts);
    blockSize = min(numMatrices, max(1, floor( ...
        options.MaxWorkingMemoryMB*2^20 / bytesPerPage)));
    estimates = zeros(numMatrices, 1);

    for first = 1:blockSize:numMatrices
        indices = first:min(first + blockSize - 1, numMatrices);
        count = numel(indices);
        startsBlock = repmat(starts, 1, 1, count);
        if options.ScreenIterations > 0
            screeningOptions = options;
            screeningOptions.MaxIterations = options.ScreenIterations;
            screened = power_norm4_multiple_starts( ...
                A(:,:,indices), normHandle, startsBlock, screeningOptions);
            screened(~isfinite(screened)) = -Inf;
            [~, ranking] = sort(screened, 1, 'descend');
            finalists = zeros(n, options.NumFinalists, count);
            for page = 1:count
                finalists(:,:,page) = starts(:,ranking(1:options.NumFinalists,page));
            end
            refined = power_norm4_multiple_starts( ...
                A(:,:,indices), normHandle, finalists, options);
            refined(~isfinite(refined)) = -Inf;
            blockEstimates = max(max(screened, [], 1), max(refined, [], 1));
        else
            candidates = power_norm4_multiple_starts( ...
                A(:,:,indices), normHandle, startsBlock, options);
            candidates(~isfinite(candidates)) = -Inf;
            blockEstimates = max(candidates, [], 1);
        end
        failed = find(~isfinite(blockEstimates), 1);
        if ~isempty(failed)
            error('optimizer:compute_norm4_fixed_starts:NoCandidate', ...
                'No start produced a finite feasible candidate for page %d.', ...
                indices(failed));
        end
        estimates(indices) = blockEstimates(:);
    end
end

function options = complete_options(options, numStarts)
    defaults = struct('MaxIterations', 1000, 'NormTolerance', 1e-10, ...
        'StationarityTolerance', 1e-8, 'FeasibilityTolerance', 1e-12, ...
        'MaxWorkingMemoryMB', 128, 'ScreenIterations', 0, 'NumFinalists', 3);
    if nargin < 1 || isempty(options)
        options = struct();
    end
    if ~isstruct(options) || ~isscalar(options)
        error('optimizer:compute_norm4_fixed_starts:InvalidOptions', ...
            'options must be a scalar struct.');
    end
    names = fieldnames(options);
    for index = 1:numel(names)
        name = names{index};
        if ~isfield(defaults, name)
            error('optimizer:compute_norm4_fixed_starts:UnknownOption', ...
                'Unknown options field: %s.', name);
        end
        defaults.(name) = options.(name);
    end
    options = defaults;
    if ~valid_integer(options.MaxIterations, 1) || ...
            ~valid_integer(options.ScreenIterations, 0) || ...
            ~valid_integer(options.NumFinalists, 1) || ...
            options.ScreenIterations > options.MaxIterations || ...
            (options.ScreenIterations > 0 && options.NumFinalists > numStarts)
        error('optimizer:compute_norm4_fixed_starts:InvalidCount', ...
            'Invalid iteration count or screening settings.');
    end
    toleranceNames = {'NormTolerance', 'StationarityTolerance', 'FeasibilityTolerance'};
    for index = 1:numel(toleranceNames)
        value = options.(toleranceNames{index});
        if ~isa(value, 'double') || ~isscalar(value) || ~isreal(value) || ...
                ~isfinite(value) || value <= 0
            error('optimizer:compute_norm4_fixed_starts:InvalidTolerance', ...
                '%s must be a finite positive double scalar.', toleranceNames{index});
        end
    end
    value = options.MaxWorkingMemoryMB;
    if ~isa(value, 'double') || ~isscalar(value) || ~isreal(value) || ...
            ~isfinite(value) || value <= 0
        error('optimizer:compute_norm4_fixed_starts:InvalidMemoryBudget', ...
            'MaxWorkingMemoryMB must be a finite positive double scalar.');
    end
end

function valid = valid_integer(value, minimum)
    valid = isa(value, 'double') && isscalar(value) && isreal(value) && ...
        isfinite(value) && value >= minimum && value == floor(value);
end
