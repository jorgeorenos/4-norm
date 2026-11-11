function banks = create_norm4_start_banks(dimension, numStarts, seed)
%CREATE_NORM4_START_BANKS Create reproducible fixed Gaussian start banks.
%   BANKS = CREATE_NORM4_START_BANKS(DIMENSION, NUMSTARTS, SEED) returns
%   independent DIMENSION-by-NUMSTARTS arrays BANKS.forward and BANKS.inverse.
%   A local stream is used, so the caller's global RNG is not read or changed.
    arguments
        dimension (1,1) double {mustBeReal, mustBeFinite, mustBeInteger, mustBePositive}
        numStarts (1,1) double {mustBeReal, mustBeFinite, mustBeInteger, mustBePositive}
        seed (1,1) double {mustBeReal, mustBeFinite, mustBeInteger, mustBeNonnegative}
    end

    stream = RandStream('mt19937ar', 'Seed', seed);
    banks = struct();
    banks.forward = nonzero_gaussian_columns(stream, dimension, numStarts);
    banks.inverse = nonzero_gaussian_columns(stream, dimension, numStarts);
    banks.seed = seed;
end

function starts = nonzero_gaussian_columns(stream, dimension, numStarts)
    starts = randn(stream, dimension, numStarts);
    zeroColumns = ~any(starts, 1);
    while any(zeroColumns)
        starts(:,zeroColumns) = randn(stream, dimension, sum(zeroColumns));
        zeroColumns = ~any(starts, 1);
    end
end
