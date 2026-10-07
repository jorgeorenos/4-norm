function angles = so3_zyx_euler_angles(rotations)
%SO3_ZYX_EULER_ANGLES Extract ZYX Euler angles from SO(3) page arrays.
%   ANGLES = SO3_ZYX_EULER_ANGLES(ROTATIONS) returns a 3-by-M array with
%   rows [alpha; beta; theta] for Q = Rz(theta)*Ry(beta)*Rx(alpha).
%   The selected identity-centered representatives avoid the gimbal-lock
%   branch for this visualization.
    if ~isa(rotations, 'double') || ~isreal(rotations) || ...
            ndims(rotations) > 3 || size(rotations,1) ~= 3 || ...
            size(rotations,2) ~= 3 || isempty(rotations) || ...
            any(~isfinite(rotations(:)))
        error('so3_zyx_euler_angles:InvalidRotation', ...
            'rotations must be a nonempty finite real 3-by-3 page array.');
    end

    betaArgument = -reshape(rotations(3,1,:), 1, []);
    beta = asin(max(-1, min(1, betaArgument)));
    alpha = atan2(reshape(rotations(3,2,:), 1, []), ...
        reshape(rotations(3,3,:), 1, []));
    theta = atan2(reshape(rotations(2,1,:), 1, []), ...
        reshape(rotations(1,1,:), 1, []));
    angles = [alpha; beta; theta];
end
