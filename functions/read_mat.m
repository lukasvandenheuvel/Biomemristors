function data = read_mat(folder, condition_string)
    % Read data that was manually checked and saved in a .mat file.
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
    
    % Initialize data structure with the correct fieldnames
    file1_data  = load(files_to_load{1});
    fields      = fieldnames(file1_data.data_clean);
    c           = cell(length(fields),1);
    data        = cell2struct(c,fields);
    % Append data from files to data
    for i=1:length(files_to_load)
        file_data = load(files_to_load{i},'data_clean');
        data = [data,file_data.data_clean];
    end
    data = data(2:end); % Remove the 1st (empty) entry
end