function value = angular_amplitude(angle, Q)
%ANGULAR_AMPLITUDE Evaluate the 4-norm amplification in a planar direction.
    direction = [cos(angle); sin(angle)];
    x = direction / norm4_columns(direction);
    value = norm4_columns(Q*x);
end
