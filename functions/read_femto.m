function data = read_femto(filepath,data)
    % This function reads a file saved by the femto program
    % --------------------------------------------------------------------
    % Check if file exists
    assert(isfile(filepath), ['The specified file ',filepath,' does not exist']);
    % Load data
    data_array = load(filepath);
    if isempty(data_array)
        return
    end
    data.V = 1e3*data_array(:,2); % conversion from V to mV
    data.I = 1e9*data_array(:,3); % conversion from A to nA
end