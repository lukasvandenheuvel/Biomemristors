function down_step_indeces = find_down_steps(V_signal, min_step_size, min_step_length, n_relax)
    assert(size(V_signal,2) == 1, 'V_signal must be a column vector!');
    % Find downward voltage jumps
    down_steps  = find(V_signal(2:end)-V_signal(1:end-1) < -min_step_size ...
                      & V_signal(2:end) < 0);
    % Find downward jumps
    up_steps    = find(V_signal(2:end)-V_signal(1:end-1) > min_step_size);

    % Match downward jumps with upward jumps
    D = up_steps' - down_steps;
    D(D<0) = inf;
    [~,min_up_inds]   = min(D);
    [~,min_down_inds] = min(D');

    min_col_array = (1:size(D,1))'==min_up_inds;
    min_row_array = min_down_inds'==(1:size(D,2));

    A = zeros(size(D));
    A(min_col_array) = A(min_col_array) + 1;
    A(min_row_array) = A(min_row_array) + 1;

    [matched_down,matched_up] = find(A==2);
    down_step_indeces = [down_steps(matched_down),up_steps(matched_up)];
    % Check if step is long enough
    if ~isempty(down_step_indeces)
        down_step_indeces(down_step_indeces(:,2)-down_step_indeces(:,1) < min_step_length,:) = [];
    end
    % Check if voltage is constant between up and down step by checking its
    % derivative halfway through
    constant = true(size(down_step_indeces,1),1);
    for s=1:size(down_step_indeces,1)
        dV = diff(V_signal(down_step_indeces(s,1)+n_relax:down_step_indeces(s,2)-n_relax));
        constant(s) = sum(abs(dV))/length(dV) < 1e-8;
    end
    down_step_indeces = down_step_indeces(constant,:);
end