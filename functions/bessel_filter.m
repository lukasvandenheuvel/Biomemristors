function x_filtered = bessel_filter(x,N,Fc,Fs)
    [z,p,k]     = besself(N,Fc);          % Bessel analog filter design
    [num,den]   = zp2tf(z,p,k);             % Convert to transfer function form
    [numd,dend] = bilinear(num,den,Fs);   % Analog to Digital conversion

    x_filtered  = filtfilt(numd,dend,x);
end