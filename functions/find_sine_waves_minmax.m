function sine_wave_indeces = find_sine_waves_minmax(V_signal, fft_peak_theshold)
    % This function returns an Nx2 matrix, where N is the number of
    % identified waves.
    % The first column are the starting indeces of an identified wave.
    % The second column are the ending indeces of an indentified wave.
    % --------------------------------------------------------------------
    % A wave is 'valid' if:
    % - it starts at a zero crossing, increases to + voltages, then
    %   decreases to - voltages, and goes back to 0.
    % - it is a sine wave.
    % --------------------------------------------------------------------
    % fft_peak_theshold   = 0.89;  % This threshold serves to identify whether a wave is a sine, using the relative height of the largest peak in the FFT. It was empirially estimated. It can't be too large since the higher-frequency recordings are sometimes not a perfect sinewave, but we still want to include them. It can't be too low because we don't want corrupted waves in our analysis. I found 0.89 to be a perfect threshold and I'm confident about it as long as we don't change the sampling frequency.
    % Find zero crossings
    %zero_crossings = ((sign(V_signal(1:end-1)) .* sign(V_signal(2:end))) == -1) ...  % Find points where the sign of the voltage flips,
    %                 | ((sign(V_signal(1:end-1)) .* sign(V_signal(2:end))) == 0);    % or where the sign is zero 
    
    min_sep = 0.02e5;
    %local_max = islocalmax(V_signal,'MinProminence',min_prominence);
    local_min = islocalmin(V_signal,'MinSeparation',min_sep);
    
    % Get individual windows (every signal in between 2 subsequent minima is a window)
    local_minima_indeces = find(local_min);
    V_windows = cell(1,length(local_minima_indeces)-1);
    for i = 1:length(local_minima_indeces)-1
        start_index  = local_minima_indeces(i);
        end_index    = local_minima_indeces(i+1);
        V_windows{i} = V_signal(start_index:end_index);
    end

    % Get individual windows (every signal in between 2 crossings is a window)
    %zero_crossing_indeces = find(zero_crossings);
    %V_windows = cell(1,length(zero_crossing_indeces)-2);
    %for i = 1:length(zero_crossing_indeces)-2
    %    start_index  = zero_crossing_indeces(i);
    %    end_index    = zero_crossing_indeces(i+2);
    %    V_windows{i} = V_signal(start_index:end_index);
    %end

    % Filter the waves that have (1) a large peak in their FFT (= they are a sine wave)
    % and (2) a negative phase at the peak frequency (= they first have a peak up, then down).
    fft_peak_score      = zeros(1,length(V_windows));
    phase               = zeros(1,length(V_windows));
    for i=1:length(V_windows)
        [FFT_amplitude, FFT_phase]  = calculate_fourier_spectrum(V_windows{i});
        fft_peak_score(i)           = max(FFT_amplitude) / sum(FFT_amplitude);      % This ratio is high if one frequency is dominant
        phase_at_max_frequencies    = FFT_phase(FFT_amplitude==max(FFT_amplitude)); % Get phase at maximum frequency/ies.
        phase(i)                    = phase_at_max_frequencies(1);                  % Take first element, just in case there is more than 1 maximal value
    end
    valid_window_indeces = find((fft_peak_score > fft_peak_theshold) ...
                                & (phase < 0));
    sine_wave_indeces = [local_minima_indeces(valid_window_indeces), ...
                          local_minima_indeces(valid_window_indeces+1)];
end