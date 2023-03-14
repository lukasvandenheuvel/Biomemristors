function [I_offset, num_pores_estimate] = fit_current_to_baseline(I_data, V_data, iv_curve_table, V_window)
    V_low  = V_window(1);
    V_high = V_window(2);
    % Obtain the parts of the curve where the pores are closing 
    % (=> Voltage increases in absolute value)
    [c1,~,c3,~]     = split_cycle([V_data,I_data], V_data);                 % split cycle in 4 parts: (1) V=0 to V=max, (2) V=max to V=0, (3) V=0 to V=min, (4) V=min to V=0
    cycle_closing   = [c3; c1];                                             % represents the part of the cycle where pores are closing
    V_closing       = cycle_closing(:,1);                               
    I_closing       = cycle_closing(:,2);
    % Obtain the part where voltage is low so that no gating takes place
    no_gating_inds       = (V_closing > V_low & V_closing < V_high);        % assume no gating takes place in the closing cycle between V_low and V_high
    V_cycle_no_gating    = V_closing(no_gating_inds);
    I_cycle_no_gating    = I_closing(no_gating_inds);
    I_bl_no_gating       = interp1(iv_curve_table.V_axis_mV, iv_curve_table.I_pore_nA, V_cycle_no_gating);
    b = fit_linear_model(I_cycle_no_gating, I_bl_no_gating);
    I_offset = b(1);
    num_pores_estimate = b(2);
end