function kappa = compute_cond4_pq_fixed_starts(P, Q, banks, normHandle, options)
%COMPUTE_COND4_PQ_FIXED_STARTS Estimate kappa_4(P*Q) from fixed start banks.
%   KAPPA = COMPUTE_COND4_PQ_FIXED_STARTS(P,Q,BANKS,NORMHANDLE,OPTIONS)
%   evaluates one condition-number estimate per page of Q. P is fixed and
%   nonsingular. BANKS.forward is reused for every P*Q page and BANKS.inverse
%   is reused for every inverse page, so evaluation order does not affect the
%   objective. The inverse is formed as Q'*(P\I), not as (P*Q)'.
    if nargin < 5
        options = struct();
    end
    if ~isa(P, 'double') || ~isreal(P) || ~ismatrix(P) || isempty(P) || ...
            size(P,1) ~= size(P,2) || any(~isfinite(P(:))) || rcond(P) == 0
        error('optimizer:compute_cond4_pq_fixed_starts:InvalidFactor', ...
            'P must be a finite real nonsingular double square matrix.');
    end
    n = size(P, 1);
    if ~isa(Q, 'double') || ~isreal(Q) || ndims(Q) > 3 || isempty(Q) || ...
            size(Q,1) ~= n || size(Q,2) ~= n || any(~isfinite(Q(:)))
        error('optimizer:compute_cond4_pq_fixed_starts:InvalidRotation', ...
            'Q must be a finite real double n-by-n matrix or page array.');
    end
    if ~isstruct(banks) || ~isscalar(banks) || ...
            ~isfield(banks, 'forward') || ~isfield(banks, 'inverse')
        error('optimizer:compute_cond4_pq_fixed_starts:InvalidBanks', ...
            'banks must contain forward and inverse start arrays.');
    end

    impact = pagemtimes(P, Q);
    pInverse = P \ eye(n);
    impactInverse = pagemtimes(Q, 'transpose', pInverse, 'none');
    forwardNorm = compute_norm4_fixed_starts( ...
        impact, normHandle, banks.forward, options);
    inverseNorm = compute_norm4_fixed_starts( ...
        impactInverse, normHandle, banks.inverse, options);
    kappa = forwardNorm .* inverseNorm;
end
