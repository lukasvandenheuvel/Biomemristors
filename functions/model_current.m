function I_modeled = model_current(input_time,input_V,ode_time,ode_n,g_data,g_closed)
    V_interp  = interp1(input_time,input_V,ode_time);
    g_interp  = interp1(g_data(:,1),g_data(:,3),V_interp);
    I_modeled = (ode_n.*g_interp + (1-ode_n)*g_closed) .* V_interp;
end