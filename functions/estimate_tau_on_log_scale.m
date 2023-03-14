function [A,tau] = estimate_tau_on_log_scale(y,x,y0,num_A_points_tested)
    % Assumed model:
    % y = A + (1-A)*exp(-x/tau)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    s  = sign(y0);
    y0 = abs(y0);
    y  = abs(y);
    w  = logspace(0,-5,length(x))';
    % Assume A is somewhere between 0.5*min(y) and 1.5*max(y)
    A_array  = linspace(0.5*min(y),1.5*max(y),num_A_points_tested); 
    R2_array = zeros(1,length(A_array)) + Inf;
    for i = 1:length(A_array)
        A = A_array(i);
        y_log = log(abs(y0-A)) - log(abs(y-A));
        [tau,S] = weighted_linear_regression(y_log,x,w);
        if isreal(tau) && tau > 0
            R2_array(i) = S;
        end
    end

    [maxR2,maxInd] = min(R2_array);
    A = A_array(maxInd);
    y_log = log(abs(y0-A)) - log(abs(y-A));
    A = s*A;
    [tau,~] = weighted_linear_regression(y_log,x,w);
end