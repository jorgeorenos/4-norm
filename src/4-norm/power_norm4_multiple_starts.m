function values = power_norm4_multiple_starts(Q, normHandle, starts, options)
%POWER_NORM4_MULTIPLE_STARTS Generalized power iteration over matrix pages.
%   Q is n-by-n-by-M; starts is n-by-S-by-M and is normalized once here
%   through normHandle. Returns an S-by-M array of best feasible amplitudes.
%   A trajectory stops independently; completed matrix pages are removed.
%   A trajectory with no feasible candidate returns -Inf.
    n = size(Q, 1);
    numStarts = size(starts, 2);
    numMatrices = size(Q, 3);
    vectorized = strcmp(func2str(normHandle), 'norm_4');
    X = starts ./ apply_norm4_handle(normHandle, starts, vectorized);
    Y = pagemtimes(Q, X);
    oldValues = apply_norm4_handle(normHandle, Y, vectorized);
    bestValues = oldValues;
    if vectorized
        feasibility = abs(sum(X.^4, 1).^(1/4) - 1);
    else
        feasibility = abs(apply_norm4_handle(normHandle, X, false) - 1);
    end
    bestValues(feasibility > options.FeasibilityTolerance) = -Inf;
    V = pagemtimes(Q, 'transpose', Y.^3, 'none');
    active = true(1, numStarts, numMatrices);
    pageIndices = 1:numMatrices;
    values = -Inf(numStarts, numMatrices);
    startIds = repmat(1:numStarts, 1, 1, numMatrices);

    for k = 1:options.MaxIterations
        u = sign(V).*abs(V).^(1/3);
        active = active & all(isfinite(u), 1) & any(u ~= 0, 1);
        % Inactive trajectories use a harmless direction so no invalid or
        % zero denominator is passed through a user's vector norm handle.
        u = safe_directions(u, active, n);
        X = u ./ apply_norm4_handle(normHandle, u, vectorized);
        active = active & all(isfinite(X), 1);
        X = safe_directions(X, active, n);
        Y = pagemtimes(Q, X);
        active = active & all(isfinite(Y), 1);
        Y = safe_directions(Y, active, n);
        newValues = apply_norm4_handle(normHandle, Y, vectorized);
        if vectorized
            feasibility = abs(sum(X.^4, 1).^(1/4) - 1);
        else
            feasibility = abs(apply_norm4_handle(normHandle, X, false) - 1);
        end
        improves = active & feasibility <= options.FeasibilityTolerance & ...
            newValues > bestValues;
        bestValues(improves) = newValues(improves);

        changeScale = max(1, max(abs(newValues), abs(oldValues)));
        active = active & newValues >= oldValues - 128*eps(changeScale);
        vNew = pagemtimes(Q, 'transpose', Y.^3, 'none');
        lambda = newValues.^4;
        active = active & isfinite(lambda) & all(isfinite(vNew), 1);
        relativeChange = abs(newValues - oldValues) ./ changeScale;
        eligible = active & relativeChange <= options.NormTolerance & ...
            feasibility <= options.FeasibilityTolerance;
        if any(eligible(:))
            % Cubes and residuals are needed only for candidates that can stop.
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
                    apply_norm4_handle(normHandle, vColumns(:,columns), vectorized) + ...
                    reshape(lambda(columns), 1, []).* ...
                    apply_norm4_handle(normHandle, xCube, vectorized), ...
                    realmin);
                relativeResidual = apply_norm4_handle( ...
                    normHandle, residual, vectorized) ./ denominator;
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
        % Discard inactive slots when all remaining pages have substantially
        % fewer live starts. Preserve original start IDs for the output.
        slots = size(active, 2);
        if vectorized && slots > 1
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
