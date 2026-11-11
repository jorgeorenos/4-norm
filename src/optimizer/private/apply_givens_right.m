function rotations = apply_givens_right(rotations, firstIndex, secondIndex, angles)
%APPLY_GIVENS_RIGHT Apply page-wise right Givens rotations.
    pageCount = size(rotations, 3);
    angles = reshape(angles, 1, 1, []);
    if isscalar(angles) && pageCount > 1
        angles = repmat(angles, 1, 1, pageCount);
    end
    if numel(angles) ~= pageCount
        error('optimizer:apply_givens_right:InvalidAngles', ...
            'Provide one angle per rotation page.');
    end
    cosine = cos(angles);
    sine = sin(angles);
    firstColumns = rotations(:,firstIndex,:);
    secondColumns = rotations(:,secondIndex,:);
    rotations(:,firstIndex,:) = firstColumns.*cosine + secondColumns.*sine;
    rotations(:,secondIndex,:) = -firstColumns.*sine + secondColumns.*cosine;
end
