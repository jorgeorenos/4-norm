function kappa = compute_4_cond(Q, normHandle, options)
%COMPUTE_4_COND Estimate the induced 4-norm condition number per matrix.
%   Q is a real n-by-n matrix or n-by-n-by-M array whose pages are assumed
%   to belong to SO(n). normHandle must implement the vector 4-norm. The
%   inverse is evaluated as Q' and no global maximum is certified. For a
%   matrix input the result is scalar; otherwise it is an M-by-1 vector.
    if nargin < 3
        options = struct();
    end
    if ~isa(normHandle, 'function_handle')
        error('compute_4_cond:InvalidHandle', ...
            'normHandle must be a function handle implementing the vector 4-norm.');
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
    values = compute_4_norm(batch, normHandle, options);
    kappa = values(1:2:end) .* values(2:2:end);
end
