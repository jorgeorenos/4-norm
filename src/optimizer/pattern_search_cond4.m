function result = pattern_search_cond4(P, initialRotations, banks, normHandle, normOptions, searchOptions)
%PATTERN_SEARCH_COND4 Refine rotations to reduce estimated kappa_4(P*Q).
%   RESULT = PATTERN_SEARCH_COND4(P,INITIALROTATIONS,BANKS,NORMHANDLE,
%   NORMOPTIONS,SEARCHOPTIONS) runs a deterministic opportunistic pattern
%   search from every page of INITIALROTATIONS. Poll directions are right
%   Givens rotations in all coordinate planes with positive and negative mesh
%   angles, so every evaluated candidate remains in SO(n).
%
%   The same BANKS.forward and BANKS.inverse starts are used for every
%   candidate. Initial and returned rotations are canonicalized toward the
%   identity. SEARCHOPTIONS.StageIterationLimits and StageRetainedCounts can
%   define a progressive schedule. The first stage uses every initial point;
%   later stages retain the lowest current estimates without resetting their
%   rotations, meshes, counters, or direction order. Pruned trajectories have
%   termination reason 'StagePruned'. No global minimum is certified. The
% function consumes no RNG. RESULT.functionCounts includes each initial and
% poll evaluation, but not the final canonical reevaluation. Set
% SEARCHOPTIONS.TrackHistory to true to retain the accepted iterate, mesh,
% and activity state after every outer iteration for visualization.
    if nargin < 5
        normOptions = struct();
    end
    if nargin < 6
        searchOptions = struct();
    end
    options = complete_search_options(searchOptions);
    [dimension, columnCount, numStarts] = size(initialRotations);
    if dimension ~= columnCount || dimension < 2 || ...
            ~isa(initialRotations, 'double') || ~isreal(initialRotations) || ...
            any(~isfinite(initialRotations(:)))
        error('optimizer:pattern_search_cond4:InvalidRotations', ...
            'initialRotations must be a finite real double n-by-n-by-M array with n >= 2.');
    end
    if ~isa(P, 'double') || ~isequal(size(P), [dimension, dimension]) || ...
            ~isreal(P) || any(~isfinite(P(:))) || rcond(P) == 0
        error('optimizer:pattern_search_cond4:InvalidFactor', ...
            'P must be a finite real nonsingular double matrix matching the rotations.');
    end
    validate_rotations(initialRotations);
    validate_stage_options(options, numStarts);

    currentRotations = canonicalize_rotations_by_trace_batch(initialRotations);
    initialValues = compute_cond4_pq_fixed_starts( ...
        P, currentRotations, banks, normHandle, normOptions);
    currentValues = initialValues;
    meshSizes = repmat(options.InitialMeshSize, numStarts, 1);
    iterations = zeros(numStarts, 1);
    functionCounts = ones(numStarts, 1);
    terminationReasons = repmat({''}, numStarts, 1);
    terminated = false(numStarts, 1);
    if options.TrackHistory
        history = struct();
        history.rotations = zeros(dimension, dimension, numStarts, ...
            options.MaxIterations + 1);
        history.values = NaN(numStarts, options.MaxIterations + 1);
        history.meshSizes = NaN(numStarts, options.MaxIterations + 1);
        history.activeAtStart = false(numStarts, options.MaxIterations + 1);
        history.accepted = false(numStarts, options.MaxIterations + 1);
        history.terminated = false(numStarts, options.MaxIterations + 1);
        history.rotations(:,:,:,1) = currentRotations;
        history.values(:,1) = currentValues;
        history.meshSizes(:,1) = meshSizes;
        history.activeAtStart(:,1) = true;
    else
        history = struct();
    end

    [firstPlanes, secondPlanes] = find(triu(true(dimension), 1));
    planeCount = numel(firstPlanes);
    directionPlanes = repelem((1:planeCount)', 2);
    directionSigns = repmat([1; -1], planeCount, 1);
    directionCount = numel(directionPlanes);
    outerIteration = 0;
    stageIndex = 1;

    while any(~terminated)
        outerIteration = outerIteration + 1;
        active = ~terminated;
        successful = false(numStarts, 1);
        shift = mod(outerIteration - 1, directionCount);
        directionOrder = circshift((1:directionCount)', -shift);

        for orderIndex = 1:directionCount
            directionIndex = directionOrder(orderIndex);
            eligible = active & ~successful & ...
                functionCounts < options.MaxFunctionEvaluations;
            if ~any(eligible)
                break
            end
            pageIndices = find(eligible);
            planeIndex = directionPlanes(directionIndex);
            angles = directionSigns(directionIndex) .* meshSizes(pageIndices);
            candidates = apply_givens_right(currentRotations(:,:,pageIndices), ...
                firstPlanes(planeIndex), secondPlanes(planeIndex), angles);
            candidateValues = compute_cond4_pq_fixed_starts( ...
                P, candidates, banks, normHandle, normOptions);
            functionCounts(pageIndices) = functionCounts(pageIndices) + 1;
            decrease = options.FunctionTolerance .* ...
                max(1, abs(currentValues(pageIndices)));
            improved = candidateValues < currentValues(pageIndices) - decrease;
            improvedIndices = pageIndices(improved);
            if ~isempty(improvedIndices)
                currentRotations(:,:,improvedIndices) = candidates(:,:,improved);
                currentValues(improvedIndices) = candidateValues(improved);
                successful(improvedIndices) = true;
            end
        end

        iterations(active) = iterations(active) + 1;
        succeeded = active & successful;
        failed = active & ~successful;
        meshSizes(succeeded) = min(options.MaxMeshSize, ...
            meshSizes(succeeded).*options.MeshExpansionFactor);
        meshSizes(failed) = meshSizes(failed).*options.MeshContractionFactor;

        reachedMesh = active & meshSizes < options.MeshTolerance;
        reachedEvaluations = active & ...
            functionCounts >= options.MaxFunctionEvaluations & ~reachedMesh;
        reachedIterations = active & ...
            iterations >= options.MaxIterations & ~reachedMesh & ~reachedEvaluations;
        terminationReasons(reachedMesh) = {'MeshTolerance'};
        terminationReasons(reachedEvaluations) = {'MaxFunctionEvaluations'};
        terminationReasons(reachedIterations) = {'MaxIterations'};
        terminated = terminated | reachedMesh | reachedEvaluations | reachedIterations;

        if ~isempty(options.StageIterationLimits) && ...
                stageIndex < numel(options.StageIterationLimits) && ...
                outerIteration >= options.StageIterationLimits(stageIndex)
            stageIndex = stageIndex + 1;
            activeIndices = find(~terminated);
            retainCount = min(options.StageRetainedCounts(stageIndex), ...
                numel(activeIndices));
            if retainCount < numel(activeIndices)
                [~, ranking] = sortrows( ...
                    [currentValues(activeIndices), activeIndices], [1, 2]);
                prunedIndices = activeIndices(ranking(retainCount + 1:end));
                terminationReasons(prunedIndices) = {'StagePruned'};
                terminated(prunedIndices) = true;
            end
        end

        if options.TrackHistory
            historyIndex = outerIteration + 1;
            history.rotations(:,:,:,historyIndex) = currentRotations;
            history.values(:,historyIndex) = currentValues;
            history.meshSizes(:,historyIndex) = meshSizes;
            history.activeAtStart(:,historyIndex) = active;
            history.accepted(:,historyIndex) = successful;
            history.terminated(:,historyIndex) = terminated;
        end

        if strcmp(options.Display, 'iter')
            fprintf('Pattern Search iteration %d: active %d, best %.12f\n', ...
                outerIteration, sum(~terminated), min(currentValues));
        end
    end

    rawFinalRotations = currentRotations;
    rawFinalValues = currentValues;
    finalRotations = canonicalize_rotations_by_trace_batch(rawFinalRotations);
    finalValues = compute_cond4_pq_fixed_starts( ...
        P, finalRotations, banks, normHandle, normOptions);

    % Finite-start estimates are not exactly invariant under signed
    % permutations. Never return a canonical estimate worse than its
    % canonical initial point under the optimization banks.
    useInitial = finalValues > initialValues;
    for index = find(useInitial).'
        finalRotations(:,:,index) = canonicalize_rotations_by_trace_batch( ...
            initialRotations(:,:,index));
    end
    finalValues(useInitial) = initialValues(useInitial);

    [bestValue, bestStartIndex] = min(finalValues);
    result = struct();
    result.bestRotation = finalRotations(:,:,bestStartIndex);
    result.bestValue = bestValue;
    result.bestStartIndex = bestStartIndex;
    result.initialRotations = canonicalize_rotations_by_trace_batch(initialRotations);
    result.initialValues = initialValues;
    result.rawFinalRotations = rawFinalRotations;
    result.rawFinalValues = rawFinalValues;
    result.finalRotations = finalRotations;
    result.finalValues = finalValues;
    result.iterations = iterations;
    result.functionCounts = functionCounts;
    result.finalMeshSizes = meshSizes;
    result.terminationReasons = terminationReasons;
    result.options = options;
    if options.TrackHistory
        historyCount = outerIteration + 1;
        history.rotations = history.rotations(:,:,:,1:historyCount);
        history.values = history.values(:,1:historyCount);
        history.meshSizes = history.meshSizes(:,1:historyCount);
        history.activeAtStart = history.activeAtStart(:,1:historyCount);
        history.accepted = history.accepted(:,1:historyCount);
        history.terminated = history.terminated(:,1:historyCount);
        history.outerIterations = outerIteration;
    end
    result.history = history;

    if strcmp(options.Display, 'final') || strcmp(options.Display, 'iter')
        fprintf(['Pattern Search completed %d starts. Initial best %.12f, ', ...
            'final best %.12f.\n'], numStarts, min(initialValues), bestValue);
    end
end

function options = complete_search_options(options)
    defaults = struct('InitialMeshSize', 0.05, 'MeshTolerance', 1e-3, ...
        'MaxMeshSize', 0.25, 'MeshExpansionFactor', 2, ...
        'MeshContractionFactor', 0.5, 'FunctionTolerance', 1e-8, ...
        'MaxIterations', 50, 'MaxFunctionEvaluations', 2000, ...
        'Display', 'final', 'TrackHistory', false, 'StageIterationLimits', [], ...
        'StageRetainedCounts', []);
    if isempty(options)
        options = struct();
    end
    if ~isstruct(options) || ~isscalar(options)
        error('optimizer:pattern_search_cond4:InvalidOptions', ...
            'searchOptions must be a scalar struct.');
    end
    names = fieldnames(options);
    for index = 1:numel(names)
        name = names{index};
        if ~isfield(defaults, name)
            error('optimizer:pattern_search_cond4:UnknownOption', ...
                'Unknown search option: %s.', name);
        end
        defaults.(name) = options.(name);
    end
    options = defaults;

    positiveNames = {'InitialMeshSize', 'MeshTolerance', 'MaxMeshSize'};
    for index = 1:numel(positiveNames)
        value = options.(positiveNames{index});
        if ~valid_scalar(value) || value <= 0
            error('optimizer:pattern_search_cond4:InvalidMesh', ...
                '%s must be a finite positive double scalar.', positiveNames{index});
        end
    end
    if options.MaxMeshSize < options.InitialMeshSize || ...
            ~valid_scalar(options.MeshExpansionFactor) || ...
            options.MeshExpansionFactor <= 1 || ...
            ~valid_scalar(options.MeshContractionFactor) || ...
            options.MeshContractionFactor <= 0 || options.MeshContractionFactor >= 1
        error('optimizer:pattern_search_cond4:InvalidMesh', ...
            'Invalid mesh size, expansion factor, or contraction factor.');
    end
    if ~valid_scalar(options.FunctionTolerance) || options.FunctionTolerance < 0
        error('optimizer:pattern_search_cond4:InvalidTolerance', ...
            'FunctionTolerance must be a finite nonnegative double scalar.');
    end
    if ~valid_integer(options.MaxIterations, 1) || ...
            ~valid_integer(options.MaxFunctionEvaluations, 1)
        error('optimizer:pattern_search_cond4:InvalidCount', ...
            'Iteration and function-evaluation limits must be positive integers.');
    end
    if ~(ischar(options.Display) && isrow(options.Display) && ...
            any(strcmp(options.Display, {'off', 'iter', 'final'})))
        error('optimizer:pattern_search_cond4:InvalidDisplay', ...
            'Display must be ''off'', ''iter'', or ''final''.');
    end
    if ~islogical(options.TrackHistory) || ~isscalar(options.TrackHistory)
        error('optimizer:pattern_search_cond4:InvalidTrackHistory', ...
            'TrackHistory must be a logical scalar.');
    end
end

function validate_stage_options(options, numStarts)
    limits = options.StageIterationLimits;
    counts = options.StageRetainedCounts;
    if isempty(limits) && isempty(counts)
        return
    end
    if isempty(limits) || isempty(counts) || ...
            ~isa(limits, 'double') || ~isrow(limits) || ~isreal(limits) || ...
            any(~isfinite(limits)) || any(limits < 1) || ...
            any(limits ~= floor(limits)) || any(diff(limits) <= 0) || ...
            limits(end) ~= options.MaxIterations || ...
            ~isa(counts, 'double') || ~isrow(counts) || ~isreal(counts) || ...
            any(~isfinite(counts)) || any(counts < 1) || ...
            any(counts ~= floor(counts)) || numel(counts) ~= numel(limits) || ...
            counts(1) ~= numStarts || any(diff(counts) >= 0)
        error('optimizer:pattern_search_cond4:InvalidStageSchedule', ...
            ['StageIterationLimits and StageRetainedCounts must be equal-length ', ...
            'double row vectors. Limits must increase to MaxIterations, and ', ...
            'counts must decrease from the number of initial rotations.']);
    end
end

function validate_rotations(rotations)
    dimension = size(rotations, 1);
    tolerance = 1e-10*max(1, dimension);
    for page = 1:size(rotations, 3)
        rotation = rotations(:,:,page);
        if norm(rotation'*rotation - eye(dimension), 'fro') > tolerance || ...
                abs(det(rotation) - 1) > tolerance
            error('optimizer:pattern_search_cond4:NotSpecialOrthogonal', ...
                'Every initial rotation must belong to SO(n) within tolerance.');
        end
    end
end

function valid = valid_scalar(value)
    valid = isa(value, 'double') && isscalar(value) && isreal(value) && isfinite(value);
end

function valid = valid_integer(value, minimum)
    valid = valid_scalar(value) && value >= minimum && value == floor(value);
end
