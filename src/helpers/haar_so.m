function Q = haar_so(n)
%HAAR_SO Generate a Haar-random matrix in the special orthogonal group.
%   Q = HAAR_SO(N) returns an N-by-N matrix in SO(N). QR factorization
%   of a Gaussian matrix, with positive diagonal in R, yields a Haar
%   sample in O(N); the final column adjustment maps it into SO(N).

    if ~isnumeric(n) || ~isscalar(n) || ~isreal(n) || ~isfinite(n) || ...
            n <= 0 || n ~= floor(n)
        error('haar_so:InvalidDimension', ...
            'n must be a finite positive integer scalar.');
    end

    if n == 1
        Q = 1;
        return;
    end

    diagonalHasZero = true;
    while diagonalHasZero
        G = randn(n, n);
        [Q, R] = qr(G);
        d = sign(diag(R));
        diagonalHasZero = any(d == 0);
    end
    Q = Q * diag(d);

    if det(Q) < 0
        Q(:, end) = -Q(:, end);
    end
end
