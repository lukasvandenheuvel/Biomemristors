function [C,n_estimate] = fit_exponential_decay_on_log_scale(I_array, time, voltage, n0, g_open, g_closed, num_A_points_tested)
    % Model:
    % ------
    % n = A + (1-A)*exp(-B*t)
    % ------
    % Take the log on both sides:
    % ln(n-A) = ln(1-A) - B*t
    % => y = ln(n-A) - ln(1-A) = b*t, where b=-B.
    % The time constant: tau = 1/B.
    % ------
    % So to find A and B, we try many different values of A, and see at
    % which value [ln(n-A) - ln(1-A)] is most linear in t (i.e. has highest
    % R2).
    % There we find b, and B=-b.
    % ------
    % Output: C = [A,tau]
    
    % Threshold to be applied to y:
    % y-values smaller than this number will skew the linear fit a lot.
    thres = 10^(-1.5);
    
    % Estimate number of pores from I at t=0.
    %[~,max_index] = max(abs(I_array));
    if I_array(end) < I_array(1)
        [~,I0_ind] = max(I_array);
    else
        [~,I0_ind] = min(I_array);
    end
    N_pores = I_array(I0_ind) / (voltage*(n0*g_open + (1-n0)*g_closed));

    % Estimate number of open pores from data
    n_estimate = (I_array/(N_pores*voltage) - g_closed) / (g_open - g_closed);

    % Grid-search on A: subtract different A-values until y is linear on a
    % y-scale.
    % So we test a number of A values, and pick the one with the highest
    % R2 on a log scale.
%     A_array = linspace(0,1,1000); % A is some value between 0 and 1
%     R2_array = zeros(1,length(A_array));
%     for i = 1:length(A_array)
%         A = A_array(i);
%         
%         y = log(abs(n_estimate-A)) - log(abs(n0-A));
%         % Remove small values as they will skew the linear fit
%         small_values = (y<log10(thres));
%         y = y(~small_values);
%         t = time(~small_values);
%         
%         [b,R2] = linearRegression(t,y);
%         if isreal(b) && b < 0
%             R2_array(i) = R2;
%         end
%     end
% 
%     [maxR2,maxInd] = max(R2_array);
%     A = A_array(maxInd);
%     y = log(abs(n_estimate-A)) - log(abs(n0-A));
%     small_values = (y<log10(thres));
%     y = y(~small_values);
%     t = time(~small_values);
%     [b,~] = linearRegression(t,y);
%     B = -b;

     [A,tau] = estimateTauOnLogScale(n_estimate,time,n0,num_A_points_tested);
    
    %{
    % Selection based on R2 and B
    if (maxR2 < 0 || B < 0) % Not an exponential at all
        % We set A to 1 and b to very small (straight line)
        A = 1;
        B = 0.3;  % NO!! We have a fitting parameter here -- talk to Simon about it
    end
    %}
    %C = [A,1/B];
    C = [A,tau];
    B = 1/tau;
    disp(['V=',num2str(voltage),', A=',num2str(A),', B=',num2str(B)]) %', R2=',num2str(maxR2)])
end