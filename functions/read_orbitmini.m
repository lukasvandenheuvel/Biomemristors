function data = read_orbitmini(folderpath,data)
    % Reads data in a folder outputted by orbit mini.
    % The function will detect the number of channels, and concatenate all
    % recordings corresponding to a channel.
    % The output (data) is a structure array with the same length as the
    % number of channels.
    % The fiels are:
    % (1) recording: contains recording (1st column is current, 2nd is voltage)
    % (2) si: sampling period
    % (3) h: header information
    % ----------------------------------------------------
    % Check if folder exists
    assert(isfolder(folderpath), ['The folder path ',folderpath,' does not exist']);
    % Read through all the files in the folder
    split_folder = split(folderpath,filesep);
    recording_specifier = split_folder{end};
    files = dir(folderpath);
    count = 0;
    channels = {};
    for f=1:length(files)
        name = files(f).name;
        if contains(name,'abf')
            count = count + 1;
            name_split = split(name, '_');
            for s = 1:length(name_split)
                snippet = name_split{s};
                found_ch = false;
                if contains(snippet,'CH')
                    channels{count} = name_split{s};
                    found_ch = true;
                    break
                end
            end
            if ~found_ch
                error('Found a .abf file with no CH indication. Aborting.')
            end
        end
    end
    % Determine which channels there were
    unique_channels = sort(sort(unique(channels)));
    num_channels    = length(unique_channels);
    % Read files channel-by-channel
    for c = 1:num_channels
        channel = unique_channels{c};
        file_specifier = [recording_specifier, '_', channel];
        % The recordings are split in files of 5 mins each. Now combine all of
        % them together.
        split_num = 0;
        file_name = [file_specifier, '_', sprintf('%03d',split_num),'.abf'];
        data(c).V = [];
        data(c).I = [];
        while isfile(fullfile(folderpath,file_name))
            % Read ABF file and append it to recording
            [d,si,h] = abf2load(fullfile(folderpath, file_name));
            data(c).I               = [data(c).I; d(:,1)];
            data(c).V               = [data(c).V; d(:,2)];
            data(c).sampling_rate   = 1/(si*1e-6);
            data(c).header          = h;
            % Increment to next part of the recording
            split_num = split_num + 1;
            file_name = [file_specifier, '_', sprintf('%03d',split_num),'.abf'];
        end
    end
end