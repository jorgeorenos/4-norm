function value = spherical_amplitude(angles, Q, normHandle)
%SPHERICAL_AMPLITUDE Evaluate 4-norm amplification in a 3-D direction.
%   ANGLES is [azimuth; polar], with polar measured from the positive z-axis.
    if ~isa(angles, 'double') || ~isvector(angles) || numel(angles) ~= 2 || ...
            ~isreal(angles) || any(~isfinite(angles(:)))
        error('spherical_amplitude:InvalidAngles', ...
            'angles must be a finite real double vector with two elements.');
    end
    azimuth = angles(1);
    polar = angles(2);
    direction = [cos(azimuth)*sin(polar); ...
        sin(azimuth)*sin(polar); cos(polar)];
    x = direction / normHandle(direction);
    value = normHandle(Q*x);
end
