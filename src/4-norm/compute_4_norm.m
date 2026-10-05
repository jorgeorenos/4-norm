function qNorm4 = compute_4_norm(Q, options)
%COMPUTE_4_NORM Estimate one induced 4-norm per square matrix.
%   Q is a real double n-by-n matrix or n-by-n-by-M array. The only
%   output is an M-by-1 vector in the same order as Q(:,:,k).
%   Each matrix uses independent random starts; no global maximum is certified.
%   Membership in SO(n) is not verified: callers working with SO(n)
%   matrices guarantee it by construction (see haar_so). The theoretical
%   bounds 1 <= ||Q||_4 <= n^(1/4) hold only for SO(n) inputs.
%   MaxWorkingMemoryMB (default 128) controls the estimated page block size.
%   ScreenIterations > 0 is an optional approximate two-stage search: run
%   every start briefly, then fully iterate only the NumFinalists strongest
%   starts per matrix. Discarded starts may otherwise have won later.
    if nargin < 2
        options = struct();
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
    for j = 1:numel(names)
        name = names{j};
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
    tolerances = {'NormTolerance', 'StationarityTolerance', 'FeasibilityTolerance'};
    for j = 1:numel(tolerances)
        tolerance = options.(tolerances{j});
        if ~isa(tolerance, 'double') || ~isscalar(tolerance) || ...
                ~isreal(tolerance) || ~isfinite(tolerance) || tolerance <= 0
            error('compute_4_norm:InvalidTolerance', ...
                '%s must be a finite positive double scalar.', tolerances{j});
        end
    end

    if ~isa(options.MaxWorkingMemoryMB, 'double') || ...
            ~isscalar(options.MaxWorkingMemoryMB) || ...
            ~isreal(options.MaxWorkingMemoryMB) || ...
            ~isfinite(options.MaxWorkingMemoryMB) || options.MaxWorkingMemoryMB <= 0
        error('compute_4_norm:InvalidMemoryBudget', ...
            'MaxWorkingMemoryMB must be a finite positive double scalar.');
    end

    % Account for matrix pages, vectors, masks, and temporary norm arrays.
    % This is a conservative working-size estimate, not a process memory cap.
    bytesPerPage = 8*(n*n + 32*n*options.NumRandomStarts + 16*options.NumRandomStarts);
    blockSize = min(numMatrices, max(1, floor( ...
        options.MaxWorkingMemoryMB*2^20 / bytesPerPage)));
    qNorm4 = zeros(numMatrices, 1);
    for first = 1:blockSize:numMatrices
        indices = first:min(first+blockSize-1, numMatrices);
        starts = randn(n, options.NumRandomStarts, numel(indices));
        columns = reshape(starts, n, []);
        zeroColumns = ~any(columns, 1);
        while any(zeroColumns)
            columns(:,zeroColumns) = randn(n, sum(zeroColumns));
            zeroColumns = ~any(columns, 1);
        end
        starts = reshape(columns, n, options.NumRandomStarts, numel(indices));
        if options.ScreenIterations == 0
            candidates = power_norm4_multiple_starts(Q(:,:,indices), starts, options);
            candidates(~isfinite(candidates)) = -Inf;
            estimates = max(candidates, [], 1);
        else
            screeningOptions = options;
            screeningOptions.MaxIterations = options.ScreenIterations;
            screened = power_norm4_multiple_starts( ...
                Q(:,:,indices), starts, screeningOptions);
            screened(~isfinite(screened)) = -Inf;
            [~, ranking] = sort(screened, 1, 'descend');
            finalists = zeros(n, options.NumFinalists, numel(indices));
            for page = 1:numel(indices)
                finalists(:,:,page) = starts(:,ranking(1:options.NumFinalists,page),page);
            end
            refined = power_norm4_multiple_starts( ...
                Q(:,:,indices), finalists, options);
            refined(~isfinite(refined)) = -Inf;
            estimates = max(max(screened, [], 1), max(refined, [], 1));
        end
        failed = find(~isfinite(estimates), 1);
        if ~isempty(failed)
            error('compute_4_norm:NoCandidate', ...
                'No start produced a finite feasible candidate for Q(:,:,%d).', indices(failed));
        end
        qNorm4(indices) = estimates(:);
    end
end

function valid = valid_integer(value, minimum)
    valid = isa(value, 'double') && isscalar(value) && isreal(value) && ...
        isfinite(value) && value >= minimum && value == floor(value);
end
