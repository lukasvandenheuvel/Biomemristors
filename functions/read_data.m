function data = read_data(filepath, device, filename_format)
    device_list = {'femto', 'axopatch', 'orbitmini'};
    assert(any(strcmp(device_list,device)), ['The device must be one of these measurement devices: ',strjoin(device_list)]);
    % Read conditions and retreive period, if it exists
    [~, filename_no_ext, ~] = fileparts(filepath);
    conditions = read_conditions(filename_no_ext,filename_format);
    period_in_filename = NaN;
    if any(strcmp(fieldnames(conditions), 'period'))
        if isnan(str2double(conditions.period))
            period_field        = conditions.period;
            digits_in_field     = regexp(period_field,'\d*','Match');
            strings_in_field    = period_field(isstrprop(period_field,'alpha'));
            prefactor           = 1;
            if ~isempty(strings_in_field)
                if isequal(strings_in_field,'ns')
                    prefactor = 1e-9;
                elseif isequal(strings_in_field,'us')
                    prefactor = 1e-6;
                elseif isequal(strings_in_field,'ms')
                    prefactor = 1e-3;
                elseif isequal(strings_in_field,'s')
                    prefactor = 1;
                else
                    error(['Period in filename included an unexpected string ',strings_in_field])
                end
            end
            period_in_filename = str2double(digits_in_field{1}) * prefactor;
        else
            period_in_filename = str2double(conditions.period);
        end
    end
    % Initialize empty data structure
    fields      = {'V','I','sampling_rate','header','conditions','device'};
    empty_data  = cell(length(fields),1);
    data        = cell2struct(empty_data, fields);
    % Read data
    if isequal(device,'femto')
        data = read_femto(filepath,data);
    elseif isequal(device,'axopatch')
        data = read_axopatch(filepath,data);
    elseif isequal(device,'orbitmini')
        data = read_orbitmini(filepath,data);
    end
    % Loop through all entries in data. Data has more than 1 entry if there
    % are multiple channels recorded (in the case of orbitmini).
    for i = 1:length(data)
        data(i).conditions = conditions;
        data(i).device     = device;
        % Retreive 
        % Assure the sampling rate is saved
        if isempty(data(i).sampling_rate)
            assert(~isnan(period_in_filename), 'The sampling period could not be retreived from the filename nor the file header. Aborting the script. If the sampling period is written in the filename, make sure "PERIOD" is present in the filename format.');
            data(i).sampling_rate = 1 / period_in_filename;
        end
        % Assure the sampling in the header is the same as the sampling rate in
        % the filename. If not, print a warning.
        if any(strcmp(fieldnames(conditions), 'period'))
            sampling_rate_filename = 1 / period_in_filename;
            sampling_rate_header   = data(i).sampling_rate;
            if abs(sampling_rate_filename - sampling_rate_header) > 1e-16
                disp(['WARNING: the sampling rate in the header (',num2str(sampling_rate_header),' Hz) does not match the sampling rate in the filename (',num2str(sampling_rate_filename), ' Hz). Check your sampling rate!']);
            end
        end
    end
end