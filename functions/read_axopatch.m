function data = read_axopatch(filepath,data)
    % This function reads a file saved by the axopatch program
    % --------------------------------------------------------------------
    % Check if file exists
    assert(isfile(filepath), ['The specified file ',filepath,' does not exist']);    
    % Read header
    fileID_utf8 = fopen(filepath, 'r', 'n', 'UTF-8');
    bytes = fread(fileID_utf8);
    unic = native2unicode(bytes);
    header = convertCharsToStrings(unic(1:2000));
    fclose(fileID_utf8);

    % Find the position where the binary file starts
    bytestrpattern = "Binary data start at position [Byte]: ";
    bytestrposition = strfind(header,bytestrpattern);
    byte_start = str2double(extract(extractBetween(header, bytestrposition, bytestrposition+strlength(bytestrpattern)+10), digitsPattern));

    % Read all header info
    header_info = {};
    header_string = "";
    num_info = 0;
    for i = 1:byte_start
        c = unic(i);
        if c==char(8592)   % if the character is "←"
            continue
        elseif c==char(10) % if the character is "↵"
            num_info = num_info+1;
            header_info{num_info} = strtrim(header_string);
            header_string = "";
        else
            header_string = header_string + convertCharsToStrings(c);
        end
    end
    
    % Get acquisition rate
    for i = 1:length(header_info)
        info_split = split(header_info{i}, ': ');
        if isequal(info_split{1},"Acquisition Rate [Hz]")
            sampling_rate = str2double(info_split{2});
        end
    end

    % Find which channels are in the file
    channelstrpattern  = "Channels are: ";
    channelsplit1 = split(header, channelstrpattern);
    channelsplit2 = split(channelsplit1{2}, ".");
    channels = strip(split(channelsplit2{1}, ","));
    assert(any(strcmp(channels,'I_Ion')), 'One of the channels in Axopatch file must be I_Ion.');
    assert(any(strcmp(channels,'Vbias_Ion')), 'One of the channels in Axopatch file must be Vbias_Ion.');

    % Read numeric data
    fileID_double = fopen(filepath, 'r');
    fseek(fileID_double, byte_start, 'bof');
    A = fread(fileID_double,'float32','ieee-be');
    fclose(fileID_double);

    nc = length(channels);
    for c = 1:length(channels)
        if isequal(channels{c},'I_Ion')     % the channel is current
            data.I = 1e-3 * A(c:nc:end);    % conversion from pA to nA
        elseif isequal(channels{c},'Vbias_Ion') % the channel is voltage
            data.V = A(c:nc:end);           % is in mV
        end
    end
    data.header = header_info;
    data.sampling_rate = sampling_rate;
    disp('Data loaded.')
end