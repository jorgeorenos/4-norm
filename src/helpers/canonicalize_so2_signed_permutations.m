function [representatives, maxTrace] = canonicalize_so2_signed_permutations(rotations)
%CANONICALIZE_SO2_SIGNED_PERMUTATIONS Center SO(2) rotations at identity.
%   REPRESENTATIVES = CANONICALIZE_SO2_SIGNED_PERMUTATIONS(ROTATIONS)
%   replaces every 2-by-2 page Q by Q*S, where S is one of the four signed
%   permutation matrices in SO(2) that maximizes trace(Q*S). Thus the
%   canonical rotation angle belongs to [-pi/4, pi/4].
    if ~isa(rotations, 'double') || ~isreal(rotations) || ...
            ndims(rotations) > 3 || size(rotations,1) ~= 2 || ...
            size(rotations,2) ~= 2 || isempty(rotations) || ...
            any(~isfinite(rotations(:)))
        error('canonicalize_so2_signed_permutations:InvalidRotation', ...
            'rotations must be a nonempty finite real 2-by-2 page array.');
    end

    signedPermutations = cat(3, eye(2), [0, -1; 1, 0], ...
        -eye(2), [0, 1; -1, 0]);
    pageCount = size(rotations, 3);
    maxTrace = -Inf(pageCount, 1);
    bestCandidate = ones(pageCount, 1);
    for candidateIndex = 1:size(signedPermutations, 3)
        candidate = pagemtimes(rotations, signedPermutations(:,:,candidateIndex));
        candidateTrace = reshape(candidate(1,1,:) + candidate(2,2,:), [], 1);
        improves = candidateTrace > maxTrace;
        maxTrace(improves) = candidateTrace(improves);
        bestCandidate(improves) = candidateIndex;
    end

    representatives = zeros(size(rotations));
    for candidateIndex = 1:size(signedPermutations, 3)
        pages = bestCandidate == candidateIndex;
        if any(pages)
            representatives(:,:,pages) = pagemtimes(rotations(:,:,pages), ...
                signedPermutations(:,:,candidateIndex));
        end
    end
end
