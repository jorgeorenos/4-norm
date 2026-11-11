function [representatives, max_trace] = canonicalize_rotations_by_trace_batch(rotations)
    arguments
        rotations {mustBeNumeric, mustBeReal, mustBeFinite}
    end
%CANONICALIZE_ROTATIONS_BY_TRACE_BATCH Choose identity-centered rotations.
%   REPRESENTATIVES = CANONICALIZE_ROTATIONS_BY_TRACE_BATCH(ROTATIONS)
%   replaces each rotation by its identity-centered signed-permutation
%   representative. [REPRESENTATIVES, MAX_TRACE] also returns its trace.
%
%   Inputs:
%       ROTATIONS       - Finite real K-by-K-by-B array of rotations in
%           SO(K), with K >= 2. Each page is one independent rotation.
%
%   Outputs:
%       REPRESENTATIVES - K-by-K-by-B array; for each input page Q, the
%           corresponding page is Q*S for a signed permutation S in SO(K)
%           maximizing tr(Q*S).
%       MAX_TRACE       - B-by-1 vector of traces of the representatives;
%           returned only when requested.
%
%   The trace maximizer minimizes the Frobenius distance to the identity:
%
%       \|Q S - I\|_F^2 = 2K - 2\operatorname{tr}(QS).
%
%   The assignment-and-cycle-correction method avoids enumerating all
%   signed permutations. Requires matchpairs (Statistics and Machine
%   Learning Toolbox).
%
%   Additional documentation:
%       docs/condition_number/norm1/condition_number-norm1.qmd
%
%   Example:
%       [q, upper] = qr(randn(3));
%       q = q .* sign(diag(upper)).';
%       if det(q) < 0, q(:,end) = -q(:,end); end
%       [representative, maximum] = canonicalize_rotations_by_trace_batch(q);
%
%   See also INITIALIZE_CONDITION_NUMBER.

    [dimension,column_count,batch_count] = size(rotations);
    assert(dimension == column_count && dimension >= 2, ...
        'rotations must be a finite real K-by-K-by-B array with K >= 2.');
    representatives = zeros(size(rotations),'like',rotations);
    if nargout > 1, max_trace = zeros(batch_count,1); else, max_trace = []; end
    row_indices = (1:dimension)';

    for draw_index = 1:batch_count
        rotation = rotations(:,:,draw_index);
        weights = abs(rotation);
        % Negative costs and costUnmatched=1 force a complete assignment.
        assignments = matchpairs(-weights,1);
        assert(size(assignments,1) == dimension,'The assignment is incomplete.');
        permutation = zeros(dimension,1);
        permutation(assignments(:,1)) = assignments(:,2);
        linear_indices = sub2ind([dimension dimension],row_indices,permutation);
        assigned_values = rotation(linear_indices);
        signs = sign(assigned_values);
        signs(signs == 0) = 1;

        % det(rotation(:,permutation).*signs.') =
        % permutation_sign(permutation)*prod(signs) for rotation in SO(K).
        if permutation_sign(permutation)*prod(signs) < 0
            % Option 1: keep p and flip the smallest assigned entry.
            [smallest_value,flip_index] = min(abs(assigned_values));
            best_cost = 2*smallest_value;
            best_permutation = permutation;
            best_signs = signs;
            best_signs(flip_index) = -best_signs(flip_index);

            % Option 2: correct the orientation through a cycle, possibly
            % followed by a sign flip; prune using nonnegative bounds.
            if best_cost > 0
                [cycle_cost,cycle_indices] = find_correcting_cycle(rotation,weights,permutation,signs,best_cost);
                if ~isempty(cycle_indices) && cycle_cost < best_cost
                    best_permutation = permutation;
                    best_permutation(cycle_indices) = permutation(cycle_indices([2:end 1]));
                    cycle_values = rotation(sub2ind([dimension dimension],row_indices,best_permutation));
                    best_signs = sign(cycle_values);
                    best_signs(best_signs == 0) = 1;
                    if permutation_sign(best_permutation)*prod(best_signs) < 0
                        [~,flip_index] = min(abs(cycle_values));
                        best_signs(flip_index) = -best_signs(flip_index);
                    end
                end
            end
            permutation = best_permutation;
            signs = best_signs;
        end

        representatives(:,:,draw_index) = rotation(:,permutation).*signs.';
        if nargout > 1, max_trace(draw_index) = trace(representatives(:,:,draw_index)); end
    end
end

function [best_cost,best_cycle] = find_correcting_cycle(rotation,weights,permutation,signs,bound)
%FIND_CORRECTING_CYCLE Find a cheaper parity-correcting assignment cycle.
%   [BEST_COST,BEST_CYCLE] = FIND_CORRECTING_CYCLE(ROTATION,WEIGHTS,
%   PERMUTATION,SIGNS,BOUND) explores cycles in ascending reduced-cost order.
%   Inputs:
%       rotation    - K-by-K input rotation in SO(K).
%       weights     - K-by-K absolute values of rotation, dimensionless.
%       permutation - K-by-1 optimal unconstrained column assignment.
%       signs       - K-by-1 signs (+1 or -1) on its assigned entries.
%       bound       - Nonnegative scalar cost of one correcting sign flip.
%   Outputs:
%       best_cost   - Scalar correction cost found (at most bound).
%       best_cycle  - Row vector of assigned row indices to cycle, or [].
%   All costs are dimensionless. Stable row-wise sorting preserves column
%   order for equal edge costs; strict improvements retain the first cycle
%   on ties. Inputs are assumed valid; no RNG or external state is modified.
%   A permutation change decomposes into disjoint cycles. Each cycle loses
%   nonnegative weight relative to the optimal assignment. One cycle that
%   changes parity suffices to correct the overall parity.

    dimension = numel(permutation);
    assigned_weights = weights(sub2ind([dimension dimension],(1:dimension)',permutation));
    edge_cost = assigned_weights - weights(:,permutation);
    % Edge row_index -> next_row assigns column permutation(next_row)
    % to row row_index.

    % Shortest-path potentials give nonnegative reduced edge costs and
    % prune partial cycles without discarding a better solution.
    potential = zeros(dimension,1);
    for iteration = 1:dimension-1
        updated_potential = min(potential + edge_cost,[],1).';
        if all(updated_potential >= potential)
            break
        end
        potential = min(potential,updated_potential);
    end
    reduced_cost = max(0,edge_cost + potential - potential.');
    [~,edge_order] = sort(reduced_cost,2);
    permuted_rotation = rotation(:,permutation);
    edge_signs = -sign(permuted_rotation).*signs.';
    zero_edges = permuted_rotation == 0;
    best_cost = bound;
    best_cycle = [];
    visited = false(dimension,1);
    cycle_path = zeros(1,dimension);

    % An explicit depth-first stack avoids one function call per explored
    % partial cycle; row-order and strict tie handling match recursive DFS.
    edge_position = ones(dimension,1);
    path_cost = zeros(dimension,1);
    path_sign = ones(dimension,1);
    path_zero = false(dimension,1);
    for root = 1:dimension
        visited(root) = true;
        cycle_path(1) = root;
        edge_position(1) = 1;
        depth = 1;
        while depth > 0
            current = cycle_path(depth);
            edge_index = edge_position(depth);
            if edge_index > dimension
                visited(current) = false;
                depth = depth - 1;
                continue
            end
            edge_position(depth) = edge_index + 1;
            next_row = edge_order(current,edge_index);
            new_cost = path_cost(depth) + reduced_cost(current,next_row);
            if new_cost >= best_cost
                visited(current) = false;
                depth = depth - 1;
                continue
            end
            if next_row == root
                if depth >= 2
                    candidate = new_cost;
                    if ~(path_zero(depth) || zero_edges(current,next_row) || ...
                            path_sign(depth)*edge_signs(current,next_row) == 1)
                        cycle_rows = cycle_path(1:depth);
                        new_diagonal = assigned_weights;
                        next_rows = cycle_rows([2:end 1]);
                        new_diagonal(cycle_rows) = weights(sub2ind([dimension dimension],cycle_rows(:), ...
                            permutation(next_rows(:))));
                        candidate = candidate + 2*min(new_diagonal);
                    end
                    if candidate < best_cost
                        best_cost = candidate;
                        best_cycle = cycle_path(1:depth);
                    end
                end
            elseif next_row > root && ~visited(next_row) && depth < dimension
                visited(next_row) = true;
                depth = depth + 1;
                cycle_path(depth) = next_row;
                edge_position(depth) = 1;
                path_cost(depth) = new_cost;
                path_sign(depth) = path_sign(depth-1)*edge_signs(current,next_row);
                path_zero(depth) = path_zero(depth-1) || zero_edges(current,next_row);
            end
        end
    end
end

function sign_value = permutation_sign(permutation)
%PERMUTATION_SIGN Return the parity of a column-index permutation.
%   Input:
%       permutation - K-by-1 vector containing each index 1:K once.
%   Output:
%       sign_value  - +1 for even, -1 for odd permutation.
    dimension = numel(permutation);
    inversions = 0;
    for row_index = 1:dimension-1
        inversions = inversions + sum(permutation(row_index)>permutation(row_index+1:dimension));
    end
    sign_value = 1-2*mod(inversions,2);
end
