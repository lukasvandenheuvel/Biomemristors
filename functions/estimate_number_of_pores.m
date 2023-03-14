function N_pores = estimate_number_of_pores(V,V_estimate,I,g_measured)
    assert(any(g_measured(:,1) == V_estimate), 'Please choose a V_estimate at which the single-pore conductance was measured.')
    % Estimate the nr of pores, assuming all pores are open at V0.
    [~,V_ind]    = min(abs(V - V_estimate));
    conductance  = I./V;
    estimated_g  = conductance(V_ind);
    N_pores      = round(estimated_g/g_measured(g_measured(:,1)==V_estimate,3));   % estimated number of pores
end