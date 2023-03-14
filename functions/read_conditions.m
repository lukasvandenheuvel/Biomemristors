function conditions = read_conditions(name,format)
    % This function reads experimental conditions from a filename.
    % The format specifies the order and names of experimental conditions.
    name_split      = split(name, '_');
    format_split    = split(format, '_');
    assert(length(name_split)>=length(format_split), ['Not all experimental conditions are specified according to the format you entered: ',format])

    conditions = struct;
    for i=1:length(format_split)
        cond  = lower(format_split{i});
        conditions.(cond) = name_split{i};
    end
    
    % The name can contain more information than the format (e.g., a time
    % added to the timestamp. In that case, bundle all extra info in the
    % name together into an 'extra' field in conditions.
    if length(name_split) > length(format_split)
        conditions.extra = strjoin(name_split(i+1:end));
    end
end
