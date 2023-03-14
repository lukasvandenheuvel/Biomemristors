function up_step_indeces = find_up_steps(V_signal, min_step_size, min_step_length, n_relax)
    assert(size(V_signal,2) == 1, 'V_signal must be a column vector!');
    % Find upward voltage jumps
    up_steps    = find(V_signal(2:end)-V_signal(1:end-1) > min_step_size ...
                      & V_signal(2:end) > 0);
    % Find downward jumps
    down_steps  = find(V_signal(1:end-1)-V_signal(2:end) > min_step_size);

    % Match upward jumps with downward jumps
    D = down_steps' - up_steps;
    D(D<0) = inf;
    [~,min_down_inds]   = min(D);
    [~,min_up_inds]     = min(D');

    min_col_array = (1:size(D,1))'==min_down_inds;
    min_row_array = min_up_inds'==(1:size(D,2));

    A = zeros(size(D));
    A(min_col_array) = A(min_col_array) + 1;
    A(min_row_array) = A(min_row_array) + 1;

    [matched_up,matched_down] = find(A==2);
    up_step_indeces = [up_steps(matched_up),down_steps(matched_down)];
    % Check if step is long enough
    if ~isempty(up_step_indeces)
        up_step_indeces(up_step_indeces(:,2)-up_step_indeces(:,1) < min_step_length,:) = [];
    end
    % Check if voltage is constant between up and down step by checking its
    % derivative halfway through
    constant = true(size(up_step_indeces,1),1);
    for s=1:size(up_step_indeces,1)
        dV = diff(V_signal(up_step_indeces(s,1)+n_relax:up_step_indeces(s,2)-n_relax));
        constant(s) = sum(abs(dV))/length(dV) < 1e-8;
    end
    up_step_indeces = up_step_indeces(constant,:);
end