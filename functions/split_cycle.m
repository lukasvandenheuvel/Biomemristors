function [d1,d2,d3,d4] = split_cycle(data_to_split,V)
    % This function splits the data into 4 parts:
    % (1) from V=0 to V=max
    % (2) from V=max to V=0
    % (3) from V=0 to V=min
    % (4) from V=min to V=0
    assert(isequal(size(data_to_split,1),size(V,1)), 'data_to_split and V do not have the same nr of rows!')
    [~,max_ind] = max(V);
    [~,min_ind] = min(V);
    half_ind = round((min_ind+max_ind)/2);
    d1 = data_to_split(1:max_ind,:);                                                      % From V=0 to V=max
    d2 = data_to_split(max_ind+1:half_ind,:);                                             % From V=max to V=0
    d3 = data_to_split(half_ind+1:min_ind,:);                                             % From V=0 to V=min
    d4 = data_to_split(min_ind+1:end,:);                                                  % From V=min to V=0
end