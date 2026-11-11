function values = power_norm4_multiple_starts(Q, normHandle, starts, options)
%POWER_NORM4_MULTIPLE_STARTS Generalized 4-power iteration over pages.
%   Q is n-by-n-by-M and starts is n-by-S-by-M. normHandle must evaluate
%   the vector 4-norm for every column and page. Returns an S-by-M array of
%   best feasible amplitudes. A trajectory with no feasible candidate is -Inf.
    n = size(Q, 1);
    numStarts = size(starts, 2);
    numMatrices = size(Q, 3);
    X = starts ./ apply_norm4_handle(normHandle, starts);
    Y = pagemtimes(Q, X);
    oldValues = apply_norm4_handle(normHandle, Y);
    bestValues = oldValues;
    feasibility = abs(apply_norm4_handle(normHandle, X) - 1);
    bestValues(feasibility > options.FeasibilityTolerance) = -Inf;
    cube = Y.*Y.*Y;
    V = pagemtimes(Q, 'transpose', cube, 'none');
    active = true(1, numStarts, numMatrices);
    pageIndices = 1:numMatrices;
    values = -Inf(numStarts, numMatrices);
    startIds = repmat(1:numStarts, 1, 1, numMatrices);

    for iteration = 1:options.MaxIterations
        u = sign(V).*exp(log(abs(V))/3);
        active = active & all(isfinite(u), 1) & any(u ~= 0, 1);
        u = safe_directions(u, active, n);
        X = u ./ apply_norm4_handle(normHandle, u);
        active = active & all(isfinite(X), 1);
        X = safe_directions(X, active, n);
        Y = pagemtimes(Q, X);
        active = active & all(isfinite(Y), 1);
        Y = safe_directions(Y, active, n);
        newValues = apply_norm4_handle(normHandle, Y);
        feasibility = abs(apply_norm4_handle(normHandle, X) - 1);
        improves = active & feasibility <= options.FeasibilityTolerance & ...
            newValues > bestValues;
        bestValues(improves) = newValues(improves);

        changeScale = max(1, max(abs(newValues), abs(oldValues)));
        active = active & newValues >= oldValues - 128*eps(changeScale);
        cube = Y.*Y.*Y;
        vNew = pagemtimes(Q, 'transpose', cube, 'none');
        lambda = newValues.^4;
        active = active & isfinite(lambda) & all(isfinite(vNew), 1);
        relativeChange = abs(newValues - oldValues) ./ changeScale;
        eligible = active & relativeChange <= options.NormTolerance & ...
            feasibility <= options.FeasibilityTolerance;
        if any(eligible(:))
            columns = find(eligible(:));
            vColumns = reshape(vNew, n, []);
            xColumns = reshape(X, n, []);
            xCube = xColumns(:,columns).^3;
            residual = vColumns(:,columns) - xCube.*reshape(lambda(columns), 1, []);
            valid = all(isfinite(residual), 1);
            active(columns(~valid)) = false;
            columns = columns(valid);
            if ~isempty(columns)
                xCube = xCube(:,valid);
                residual = residual(:,valid);
                denominator = max( ...
                    apply_norm4_handle(normHandle, vColumns(:,columns)) + ...
                    reshape(lambda(columns), 1, []).* ...
                    apply_norm4_handle(normHandle, xCube), realmin);
                relativeResidual = apply_norm4_handle(normHandle, residual) ./ denominator;
                active(columns(relativeResidual(:) <= options.StationarityTolerance)) = false;
            end
        end
        V = vNew;
        oldValues = newValues;

        keep = reshape(any(active, 2), 1, []);
        if any(~keep)
            for page = find(~keep)
                validSlots = startIds(1,:,page) > 0;
                positions = startIds(1,validSlots,page) + ...
                    (pageIndices(page)-1)*numStarts;
                values(positions) = bestValues(1,validSlots,page);
            end
            pageIndices = pageIndices(keep);
            if isempty(pageIndices)
                return;
            end
            Q = Q(:,:,keep);
            V = V(:,:,keep);
            oldValues = oldValues(:,:,keep);
            bestValues = bestValues(:,:,keep);
            active = active(:,:,keep);
            startIds = startIds(:,:,keep);
        end

        slots = size(active, 2);
        if slots > 1
            counts = reshape(sum(active, 2), 1, []);
            packedSlots = max(counts);
            if packedSlots <= floor(0.75*slots)
                pages = numel(pageIndices);
                packedV = zeros(n, packedSlots, pages);
                packedOld = zeros(1, packedSlots, pages);
                packedBest = -Inf(1, packedSlots, pages);
                packedIds = zeros(1, packedSlots, pages);
                packedActive = false(1, packedSlots, pages);
                for page = 1:pages
                    validSlots = startIds(1,:,page) > 0;
                    positions = startIds(1,validSlots,page) + ...
                        (pageIndices(page)-1)*numStarts;
                    values(positions) = bestValues(1,validSlots,page);
                    live = find(active(1,:,page));
                    count = numel(live);
                    packedV(:,1:count,page) = V(:,live,page);
                    packedOld(1,1:count,page) = oldValues(1,live,page);
                    packedBest(1,1:count,page) = bestValues(1,live,page);
                    packedIds(1,1:count,page) = startIds(1,live,page);
                    packedActive(1,1:count,page) = true;
                end
                V = packedV;
                oldValues = packedOld;
                bestValues = packedBest;
                startIds = packedIds;
                active = packedActive;
            end
        end
    end

    for page = 1:numel(pageIndices)
        validSlots = startIds(1,:,page) > 0;
        positions = startIds(1,validSlots,page) + ...
            (pageIndices(page)-1)*numStarts;
        values(positions) = bestValues(1,validSlots,page);
    end
end

function X = safe_directions(X, active, n)
    if any(~active(:))
        originalSize = size(X);
        columns = reshape(X, n, []);
        inactive = ~active(:);
        columns(:,inactive) = 0;
        columns(1,inactive) = 1;
        X = reshape(columns, originalSize);
    end
end
