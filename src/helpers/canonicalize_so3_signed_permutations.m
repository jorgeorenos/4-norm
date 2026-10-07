function [representatives, maxTrace] = canonicalize_so3_signed_permutations(rotations)
%CANONICALIZE_SO3_SIGNED_PERMUTATIONS Center SO(3) rotations at identity.
%   REPRESENTATIVES = CANONICALIZE_SO3_SIGNED_PERMUTATIONS(ROTATIONS)
%   replaces each 3-by-3 page Q by Q*S, where S is the signed permutation
%   in SO(3) that maximizes trace(Q*S). The 24 candidates are enumerated,
%   avoiding a dependency on matchpairs for this SO(3)-specific operation.
    if ~isa(rotations, 'double') || ~isreal(rotations) || ...
            ndims(rotations) > 3 || size(rotations,1) ~= 3 || ...
            size(rotations,2) ~= 3 || isempty(rotations) || ...
            any(~isfinite(rotations(:)))
        error('canonicalize_so3_signed_permutations:InvalidRotation', ...
            'rotations must be a nonempty finite real 3-by-3 page array.');
    end

    candidatePermutations = perms(1:3);
    signTriples = [1, 1, 1; 1, 1, -1; 1, -1, 1; 1, -1, -1; ...
        -1, 1, 1; -1, 1, -1; -1, -1, 1; -1, -1, -1];
    candidateCount = 24;
    permutations = zeros(candidateCount, 3);
    signs = zeros(candidateCount, 3);
    candidateIndex = 0;
    for permutationIndex = 1:size(candidatePermutations,1)
        permutation = candidatePermutations(permutationIndex,:);
        permutationParity = permutation_sign(permutation);
        for signIndex = 1:size(signTriples,1)
            candidateSigns = signTriples(signIndex,:);
            if permutationParity*prod(candidateSigns) == 1
                candidateIndex = candidateIndex + 1;
                permutations(candidateIndex,:) = permutation;
                signs(candidateIndex,:) = candidateSigns;
            end
        end
    end

    pageCount = size(rotations, 3);
    bestCandidate = ones(pageCount, 1);
    maxTrace = -Inf(pageCount, 1);
    for candidateIndex = 1:candidateCount
        permutation = permutations(candidateIndex,:);
        candidateSigns = signs(candidateIndex,:);
        candidateTrace = candidateSigns(1)*reshape(rotations(1,permutation(1),:), [], 1) + ...
            candidateSigns(2)*reshape(rotations(2,permutation(2),:), [], 1) + ...
            candidateSigns(3)*reshape(rotations(3,permutation(3),:), [], 1);
        improves = candidateTrace > maxTrace;
        maxTrace(improves) = candidateTrace(improves);
        bestCandidate(improves) = candidateIndex;
    end

    representatives = zeros(size(rotations));
    for candidateIndex = 1:candidateCount
        pages = bestCandidate == candidateIndex;
        if ~any(pages)
            continue;
        end
        permutation = permutations(candidateIndex,:);
        candidateSigns = signs(candidateIndex,:);
        for column = 1:3
            representatives(:,column,pages) = ...
                candidateSigns(column)*rotations(:,permutation(column),pages);
        end
    end
end

function signValue = permutation_sign(permutation)
    inversions = 0;
    for index = 1:numel(permutation)-1
        inversions = inversions + sum(permutation(index) > permutation(index+1:end));
    end
    signValue = 1 - 2*mod(inversions, 2);
end
