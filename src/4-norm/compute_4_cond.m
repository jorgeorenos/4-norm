function kappa = compute_4_cond(Q, options)
%COMPUTE_4_COND Estimate the induced 4-norm condition number per matrix.
%   Q is a real double n-by-n matrix or n-by-n-by-M array whose pages are
%   assumed to belong to SO(n); membership is not verified. The inverse is
%   taken as the transpose, which is exact for SO(n). The only output is
%   an M-by-1 vector of estimates kappa(k) = ||Q(:,:,k)||_4*||Q(:,:,k)'||_4
%   in the order of Q(:,:,k); for one matrix the output is a scalar.
%   Pages are interleaved as Q, Q', Q, Q', ... so that a batched call
%   consumes the same random stream as sequential per-matrix calls.
%   Both factors are lower bounds of the true maxima, so the product
%   underestimates the condition number; no global maximum is certified.
%   For SO(n) inputs the true value satisfies 1 <= kappa <= n^(1/2).
    if nargin < 2
        options = struct();
    end
    if ~isa(Q, 'double') || ~isreal(Q) || ndims(Q) > 3 || ...
            isempty(Q) || size(Q,1) ~= size(Q,2) || any(~isfinite(Q(:)))
        error('compute_4_cond:InvalidMatrix', ...
            'Q must be a nonempty finite real double n-by-n or n-by-n-by-M array.');
    end
    numMatrices = size(Q, 3);
    batch = zeros(size(Q,1), size(Q,2), 2*numMatrices);
    batch(:,:,1:2:end) = Q;
    batch(:,:,2:2:end) = permute(Q, [2 1 3]);
    values = compute_4_norm(batch, options);
    kappa = values(1:2:end) .* values(2:2:end);
end
