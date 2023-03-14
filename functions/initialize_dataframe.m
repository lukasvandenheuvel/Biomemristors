function df = initialize_dataframe(fields)
    c  = cell(length(fields),1);
    df = cell2struct(c,fields);
end
