function values = evaluate_condition_number_batch(rotations, normHandle, options)
%EVALUATE_CONDITION_NUMBER_BATCH Evaluate condition estimates without RNG advance.
%   VALUES = EVALUATE_CONDITION_NUMBER_BATCH(ROTATIONS, NORMHANDLE, OPTIONS)
%   delegates to the handle-based COMPUTE_4_COND while restoring the global
%   random-number-generator state when it returns or errors. This keeps Haar
%   draws in INITIALIZE_CONDITION_NUMBER independent of the stochastic
%   multistart estimator.
    rngState = rng;
    restoreRng = onCleanup(@() rng(rngState)); %#ok<NASGU>
    values = compute_4_cond(rotations, normHandle, options);
end
