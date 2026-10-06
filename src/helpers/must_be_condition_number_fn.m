function must_be_condition_number_fn(conditionNumberFn)
%MUST_BE_CONDITION_NUMBER_FN Validate the initialization callback type.
%   The callback output is validated by initialize_condition_number because
%   its required length depends on the current sampling batch.
    if ~isa(conditionNumberFn, 'function_handle')
        error('condition_number:initialize_condition_number:InvalidCallback', ...
            'condition_number_fn must be a function handle.');
    end
end
