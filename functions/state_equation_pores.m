function dndt = state_equation_pores(t, n, timeV, V, r_open_mean, r_close_mean, V_amplitudes)
    %p_opening = @sigmoid;
    %p_closing = @sigmoid; 
    V = interp1(timeV, V, t); % Interpolate the data set (It, I) at times t
    r_open  = interp1(V_amplitudes,r_open_mean,V);
    r_close = interp1(V_amplitudes,r_close_mean,V);

    %p_opening = r * 0.5*(1+cos(pi * V));
    %p_closing = r * 0.5*(1-cos(pi * V));
    dndt = (1 - n) * r_open - n * r_close;
    %dndt = (1 - n) * p_opening - n * p_closing;
    %dndt = (1 - n) * r*p_opening(-abs(V),alpha,-x0) /0.2 - n * r*p_closing(abs(V),alpha,x0) /0.2;
end