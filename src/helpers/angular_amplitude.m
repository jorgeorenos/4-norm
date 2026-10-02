function value = angular_amplitude(angle, Q, normHandle)
%ANGULAR_AMPLITUDE Evaluate the 4-norm amplification in a planar direction.
    direction = [cos(angle); sin(angle)];
    x = direction / normHandle(direction);
    value = normHandle(Q*x);
end
