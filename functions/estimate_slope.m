function slope = estimate_slope(cycle, V_low)
    % Estimates the slope at low voltage.
    % cycles is an Nx2 array, with 1st column voltage and 2nd column
    % current.
    % V_low is a float, indicating the (positive) voltage where pores are
    % not gating yet.
    up_sweep        = cycle(round(1:end/4),:);
    up_sweep_low_V  = up_sweep(up_sweep(:,1) > 0 & up_sweep(:,1) < V_low,:);
    [b,~]           = linear_regression(up_sweep_low_V(:,1),up_sweep_low_V(:,2), true);
    slope           = b(2);
end