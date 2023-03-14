function [FFT_amplitude, FFT_phase] = calculate_fourier_spectrum(signal)
    % This function calculates the Fourier spectrum of a signal.
    % It returns the single-sided amplitude spectrum and the single-sided
    % phase spectrum.
    % --------------------------------------------------------------------
    L   = length(signal);
    FFT = fft(signal - mean(signal));                   % Matlab's built-in Fourier transform. Subtract the mean to get rid of the DC component.
    FFT_single_sided = FFT(1:round(L/2)+1);             % Take only right-half of spectrum (positive frequencies)
    FFT_amplitude    = abs(FFT_single_sided/L);         % Fourier amplitude = absolute value
    FFT_amplitude(2:end-1) = 2*FFT_amplitude(2:end-1);  % Add the fft of the negative frequencies to the single-sided spectrum
    FFT_phase        = angle(FFT_single_sided);         % The phase is the angle
end