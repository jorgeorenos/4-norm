function [best_rotation, best_condition_number, best_index] = initialize_condition_number(K, n_sample, batch_size, condition_number_fn, n_initial_condition)
    arguments
        K (1,1) {mustBeNumeric, mustBeReal, mustBeFinite, ...
            mustBeInteger, mustBeGreaterThanOrEqual(K, 2)}
        n_sample (1,1) {mustBeNumeric, mustBeReal, mustBeFinite, ...
            mustBeInteger, mustBePositive}
        batch_size (1,1) {mustBeNumeric, mustBeReal, mustBeFinite, ...
            mustBeInteger, mustBePositive}
        condition_number_fn {must_be_condition_number_fn}
        n_initial_condition (1,1) {mustBeNumeric, mustBeReal, ...
            mustBeFinite, mustBeInteger, mustBePositive, ...
            mustBeLessThanOrEqual(n_initial_condition, n_sample)} = 1
    end
%INITIALIZE_CONDITION_NUMBER Find sampled SO(K) starting rotations.
%   [BEST_ROTATION, BEST_CONDITION_NUMBER, BEST_INDEX] =
%   INITIALIZE_CONDITION_NUMBER(K, N_SAMPLE, BATCH_SIZE,
%   CONDITION_NUMBER_FN) scores batches of raw Haar rotations, selects the
%   sampled minimum, and then returns its identity-centered representative.
%   INITIALIZE_CONDITION_NUMBER(..., N_INITIAL_CONDITION) returns the
%   lowest-valued N_INITIAL_CONDITION draws instead of only the best one.
%
%   Inputs:
%       K                   - Integer dimension of SO(K), K >= 2.
%       N_SAMPLE            - Positive integer number of Haar draws.
%       BATCH_SIZE          - Positive integer maximum draws per batch.
%       CONDITION_NUMBER_FN - Function handle called once per sampling batch,
%           in draw order, on raw (not canonicalized) K-by-K-by-B rotations.
%           B <= BATCH_SIZE; for B = 1 the input is a K-by-K matrix. Returns
%           a B-element numeric row or column vector of finite, real,
%           nonnegative, dimensionless scores (a scalar for B = 1).
%           Scores must satisfy f(Q*S) = f(Q) for every signed permutation
%           S with det(S) = +1, and be independent of batch composition.
%           These properties are required, not checked. The callback must
%           not consume or modify the global RNG for batch-size invariance.
%           For example, @(rotations) condition_number(P, rotations, @norm_1)
%           returns the 1-norm condition numbers of P*Q for all input pages.
%       N_INITIAL_CONDITION - Optional positive integer between 1 and
%           N_SAMPLE; defaults to 1.
%
%   Outputs:
%       BEST_ROTATION         - K-by-K matrix when N_INITIAL_CONDITION = 1;
%           otherwise K-by-K-by-N_INITIAL_CONDITION array of representatives.
%       BEST_CONDITION_NUMBER - N_INITIAL_CONDITION-by-1 double vector of
%           callback scores on the raw draws, sorted ascending (a scalar
%           when N_INITIAL_CONDITION = 1). Scores are not recomputed after
%           canonicalization; roundoff or approximate norm evaluation can
%           cause small differences if evaluated on the returned matrices.
%       BEST_INDEX            - Corresponding 1-based indices into the draw
%           sequence, sorted by value and then index; a scalar by default.
%
%   Canonicalization minimizes the Frobenius distance from the identity
%   over signed permutations S with det(S) = +1, equivalently maximizing
%
%       \|Q S - I\|_F^2 = 2K - 2\operatorname{tr}(QS).
%
%   Equal callback values favor earlier draws. Only the selected rotations
%   are canonicalized, in batches no larger than BATCH_SIZE. Ranking uses
%   scores and indices; matrix pages are gathered once per accepted batch,
%   not shifted at every insertion. Storage is O(K^2*(BATCH_SIZE +
%   N_INITIAL_CONDITION)), including temporary selection arrays. This is a
%   sampled minimum, not a certificate of global optimality.
%
%   The global RNG state advances; set rng(seed) for reproducibility. An
%   invalid callback result raises
%   condition_number:initialize_condition_number:InvalidCallbackResult.
%   Requires matchpairs (Statistics and Machine Learning Toolbox).
%
%   Example:
%       P = [2 0; 0 1];
%       rng(245); [q, value, index] = initialize_condition_number( ...
%           2, 100, 20, @(rotation) condition_number(P, rotation, @norm_1));
%
%   See also CONDITION_NUMBER, CANONICALIZE_ROTATIONS_BY_TRACE_BATCH.

    best_rotation = zeros(K, K, 0);
    best_condition_number = zeros(0, 1);
    best_index = zeros(0, 1);
    selected_count = 0;

    for batch_start = 1:batch_size:n_sample
        count = min(batch_size, n_sample - batch_start + 1);
        rotations = zeros(K, K, count);
        for draw_index = 1:count
            [rotation, upper] = qr(randn(K), 0);
            diagonal_signs = sign(diag(upper));
            diagonal_signs(diagonal_signs == 0) = 1;
            rotation = rotation .* diagonal_signs.';
            if det(rotation) < 0
                rotation(:, end) = -rotation(:, end);
            end
            rotations(:, :, draw_index) = rotation;
        end
        values = condition_number_fn(rotations);
        if ~isnumeric(values) || ~isreal(values) || ~isvector(values) || ...
                numel(values) ~= count || any(~isfinite(values(:))) || any(values(:) < 0)
            error('condition_number:initialize_condition_number:InvalidCallbackResult', ...
                ['condition_number_fn must return one finite real nonnegative ', ...
                'numeric score per input page, as a row or column vector.']);
        end
        values = double(values(:));

        % MINK preserves input order for equal values. Earlier retained
        % draws precede the new batch, whose pages remain in draw order.
        next_count = min(n_initial_condition, selected_count + count);
        [next_values, positions] = mink( ...
            [best_condition_number(1:selected_count); values], next_count);
        from_batch = positions > selected_count;
        if ~any(from_batch)
            continue
        end
        next_rotation = zeros(K, K, next_count);
        next_index = zeros(next_count, 1);
        next_rotation(:, :, ~from_batch) = best_rotation(:, :, positions(~from_batch));
        next_index(~from_batch) = best_index(positions(~from_batch));
        batch_positions = positions(from_batch) - selected_count;
        next_rotation(:, :, from_batch) = rotations(:, :, batch_positions);
        next_index(from_batch) = batch_start + batch_positions - 1;
        best_rotation = next_rotation;
        best_condition_number = next_values;
        best_index = next_index;
        selected_count = next_count;
    end

    % Do not send representatives back to the callback: the invariant scores
    % above define the ranking, including floating-point ties.
    for batch_start = 1:batch_size:n_initial_condition
        batch_end = min(batch_start + batch_size - 1, n_initial_condition);
        best_rotation(:, :, batch_start:batch_end) = ...
            canonicalize_rotations_by_trace_batch(best_rotation(:, :, batch_start:batch_end));
    end
end
