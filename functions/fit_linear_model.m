function b = fit_linear_model(I_data, I_baseline)
    % Assume linear model: I_data = I_offset + N_pores * I_baseline
    % Here we fit I_data = b1 + b2*I_baseline, 
    % so b(1) = offset and b(2) = N_pores.
    assert(isequal(size(I_data),size(I_baseline)), 'I_data and I_baseline are not the same length!')
    X = [ones(length(I_baseline),1), I_baseline];
    b = X\I_data;
end