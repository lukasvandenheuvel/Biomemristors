function clipped = find_clipped_regions(signal, clipping_length_threshold)
    % This function returns a boolean array of the same length as signal.
    % It is 1 if the signal is 'clipped', i.e. nonchanging for at least 
    % clipping_length_threshold indeces.
    % --------------------------------------------------------------------
    derivative = diff(signal);

    x = (derivative == 0);                          % Label regions where the signal is nonchanging.
    f = find(diff([0,x',0]==1));                    % Find locations where one region flips into the other
    start_indeces   = f(1:2:end-1);                 % Start indices
    end_indeces     = f(2:2:end);                   % End indices
    counts          = end_indeces - start_indeces;  % Consecutive ones’ counts

    start_indeces = start_indeces(counts > clipping_length_threshold) - 5;  % Start indeces of each nonchanging region in signal. -5 to also take the sides 
    end_indeces   = end_indeces(counts > clipping_length_threshold) + 5;    % End indeces of each nonchanging region in signal. +5 to also take the sides 

    % make sure start_indeces and end_indeces contain existing indeces
    start_indeces = max(start_indeces,1);
    end_indeces   = min(end_indeces, length(signal));
    
    % Collect all clipped indeces in a large array
    clipped_indx  = [];
    for i = 1:length(start_indeces)
        clipped_indx = [clipped_indx, start_indeces(i):end_indeces(i)];
    end
    % Convert to a boolean
    clipped = boolean(zeros(size(signal)));
    clipped(clipped_indx) = true;
end