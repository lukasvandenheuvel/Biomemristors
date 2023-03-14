function [tot_t, tot_V, tot_I] = concatenate_femto_files(folder, condition_string)
    % This function reads all files corresponding to a single condition,
    % and then concatenates all data into three long signals (one for time,
    % one for voltage and one for current)
    % --------------------------------------------------------------------
    % Check if folder exists
    if ~isfolder(folder)
        error(['The specified folder ',folder,' does not exist.'])
    end
    % Collect files with the corresponding condition in a cell array
    files_to_load   = {};
    num_files       = 1;
    file_list       = dir(folder);
    for i=1:length(file_list)
        file = file_list(i);
        if contains(file.name,condition_string)
            files_to_load{num_files} = fullfile(file.folder, file.name);
            num_files = num_files + 1;
        end
    end
    
    if isempty(files_to_load)
        error(['No files starting with ',condition_string,' were found in ',folder,'. Are you sure you entered data_folder and condition_string correctly?'])
    end

    % Load the files and put the content in a large trace of time, current
    % and voltage.
    tot_t = [];
    tot_V = [];
    tot_I = [];
    current_time = 0;
    for i=1:length(files_to_load)
        data = load(files_to_load{i});
        if isempty(data)
            continue
        end
        tot_t = [tot_t; data(:,1) + current_time];
        tot_V = [tot_V; data(:,2)];
        tot_I = [tot_I; data(:,3)];
        current_time = tot_t(end);
    end
end