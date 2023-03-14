clear all
close all
addpath('./functions')

%% Intro
% Hello! With this function you will plot the sine waves, up-steps and
% down-steps which you selected with the script biomems_check_cycles.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% SET PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% input and output
input_path              = '/Users/lukasvandenheuvel/Documents/EPFL/MA3/LBM/Biomemdata_cleaned/Ael-K238AK242A';
conductivity_data_path  = 'none';
plotting_path           = '/Users/lukasvandenheuvel/Documents/EPFL/MA3/LBM/Biomemdata_plots';
filename_includes       = 'sine_f0.1Hz_C0_M1_Cm52.8_16-01-2023_17-01-06';          % All filenames with this string included will be combined
output_title            = 'filename_format';     % specify a title for output plots OR enter 'filename_format' to use the format of the first loaded datafile.

%input_path              = 'C:\Users\smayer\Desktop\Biomemdata-checked\Ael-wt';
%conductivity_data_path  = 'none';  % Enter a path to an xlsx file, or 'none'. If 'none', square waves will be fitted by an exponential, otherwise with the 2-state model (for which r_open and r_close will be estimated).
%plotting_path           = 'C:\Users\smayer\Desktop\Biomemdata-plots\newstuff';
%filename_includes       = 'Aelwt_KCl_1M_pH6p2_10C_100us__230203141857'; % All filenames with this string will be combined
%output_title            = filename_includes;     % specify a title for output plots OR enter 'filename_format' to use the format of the first loaded datafile.

% pore-specific parameters
g_open                  = 0.1e-9;   % rougly the conductance of a single pore when open (in Siemens)
g_closed                = 0.01e-9;  % the conductance of a single pore when closed (in Siemens)
V_low                   = 70;       % Voltage (in mV) below which almost no gating takes place. Is used to estimate number of pores. Matters for slope calculation.

% analysis parameters
num_cycles              = 15;       % number of cycles to include in analysis. The cycles with the highest maximal current will be picked.
min_num_pores           = 1;        % minimal (estimated) number of pores needed to take measurement into account
combine_cyles           = 'sum';    % How to combine multiple cycles. Choose between 'avg' and 'sum'.
normalize_by            = 'max';    % How to normalize current traces. Choose between 'max', 'slope' and 'none'.

% figure parameters
format                  =  'pdf';                   % figure extension png takes a long time
export_plots            =  true;                    % if true plots are saved in the specified format
figure_width            = 15; % inches
figure_height           = 10; % inches
figure_fontsize         = 11;
fill_color              = [0.9, 0.9, 0.9]; % color for errorbars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% START SCRIPT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load data
data = read_mat(input_path, filename_includes);
if isempty(data)
    disp(['No data found folder ',input_path,' corresponding to the name ',filename_includes,'!'])
    return
end

% conductivity data
g_data = 'none';
if ~isequal(conductivity_data_path,'none')
    assert(isfile(conductivity_data_path),['The specified conductivity path ',conductivity_data_path,' does not exist. If no conductivity data is available, enter "none".'])
    g_data = xlsread(conductivity_data_path);
    g_data(:,3) = g_data(:,2)./g_data(:,1);
end

%% Output filename and figure title
if isequal(output_title,'filename_format')
    fields = fieldnames(data(1).conditions);
    figure_title_format = "";
    for f = 1:length(fields)
        field = fields{f};
        figure_title_format = figure_title_format + data(1).conditions.(field) + " ";
    end
else
    figure_title_format = output_title + " ";
end
figure_title_format = figure_title_format + combine_cyles;
disp("Figure title format: '"+figure_title_format+"'")

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SINE WAVES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Combine sine_waves of all channels
sine_waves = initialize_dataframe(fieldnames(data(1).sine_waves));
sine_waves.maxI = [];
for i = 1:length(data)
    new_sine_waves = data(i).sine_waves;
    % Obtain max-I values
    for j = 1:length(new_sine_waves)
        new_sine_waves(j).maxI = max(new_sine_waves(j).I_cycle);
    end
    sine_waves = [sine_waves, new_sine_waves];
end
sine_waves = sine_waves(2:end); % Remove first empty field

%% ID vector that assigns each combination of a & f an integer
af = [[sine_waves(:).amp]', [sine_waves(:).freq]'];
[~,~,ID] = unique(af(:,:),'rows','stable');

%% Go through every unique combi of a & f
for i = 1:length(unique(ID)) % goes through every unique combi of a & f
    n = 0;                   % sets back number of cycles for next iteration
    unique_af_indeces = find(ID == i);                                      % lists all instances of a given unique combi of a & f
    sine_waves_af     = sine_waves(unique_af_indeces);
    % Filter the sine wave cycles with the highest amplitude
    if num_cycles < length(sine_waves_af)
        [~, ind] = sort([sine_waves_af(:).maxI]);
        sine_waves_af = sine_waves_af(ind <= num_cycles);
        disp(['Using the ',num2str(num_cycles), ' cycles with the highest amplitude.'])
    else
        disp(['Warning: you set num_cycles to ',num2str(num_cycles), ', but there are only ', num2str(length(sine_waves)), ' cycles available. Using all of them now.'])
    end
    % Find maximal cycle size (to be used for padding later on)
    max_cycle_size = 0;
    for j = 1:length(sine_waves_af)
        max_cycle_size = max(max_cycle_size, size(sine_waves_af(j).V_cycle,1));
    end
    amp               = sine_waves_af(1).amp;
    freq              = sine_waves_af(1).freq;
    V_cycles          = zeros(max_cycle_size, length(sine_waves_af));       % Store V of each cycle in one column
    I_cycles          = zeros(max_cycle_size, length(sine_waves_af));       % Store I of each cycle in one column
    cycles            = zeros(max_cycle_size, 2);                           % One column for voltage, one for current
    for j = 1:length(sine_waves_af)
        cycle = [sine_waves_af(j).V_cycle, sine_waves_af(j).I_cycle];

        % used for padding in case cycle arrays are not same length
        if size(cycle,1) < max_cycle_size
            pad = size(cycles,1) - size(cycle,1);
            cycle = [cycle;zeros(pad,2)];
        end
                
        % Estimate slope at first upwards sweep, for low voltage
        slope = 1e-6 * estimate_slope(cycle,V_low); % 1e-6 is for unit conversion
        
        % Do not include this cycle if there were too few pores
        if slope < min_num_pores*g_open
            disp(['Excluded cycle with g0 = ',num2str(slope), 'S. This conductance is too low.'])
            continue
        end
        
        n = n + 1;                              % indicates how many cycles are summed
        cycles = cycles(:,1:2) + cycle(:,1:2);  % sums the values of all cycles
        V_cycles(:,n) = cycle(:,1);
        if isequal(normalize_by,'max')
            I_cycles(:,n) = cycle(:,2) / max(abs(cycle(:,2)));  % normalize the current by max
            y_label = 'norm.I';
        elseif isequal(normalize_by,'slope')
            I_cycles(:,n) = cycle(:,2) / (1e9*slope);           % normalize the current by slope.
            y_label = 'norm.I';
        elseif isequal(normalize_by,'none')
            I_cycles(:,n) = cycle(:,2);                         % No normalization
            y_label = 'I (nA)';
        else
            error('Incorrect value for normalize_by. Choose between "max" and "slope"');
        end
    end
    
    V_cycles = V_cycles(:,1:n);
    I_cycles = I_cycles(:,1:n);
    av_voltage = cycles(:,1) / n; % averaging the voltage
    
    if isequal(combine_cyles,'sum')  % combine cycles by summing. This weighs cycles with higher current more!
        if isequal(normalize_by, 'max')
            norm_current    = cycles(:,2) / max(abs(cycles(:,2))); % normalizing summed current by its max
        elseif isequal(normalize_by,'slope')
            % Estimate slope at first upwards sweep, for low voltage
            slope           = 1e-6 * estimate_slope([av_voltage,cycles(:,2)],V_low);
            norm_current    = cycles(:,2) / (1e9*slope); % normalizing by slope. 
        elseif isequal(normalize_by,'none')
            norm_current    = cycles(:,2);
        else
            error('Incorrect value for normalize_by. Choose between "max" , "slope" or "none"');
        end
        I_cycles = norm_current; % sum of the current, normalized
    end
    
    if n==0
        disp("No valid waves to be plotted for A="+amp+" and f="+freq)
        continue
    end
    disp("Plotting sine waves for A="+amp+" and f="+freq)
    
    % Average data over cycles
    avg_voltage     = mean(V_cycles,2);   % mean over individual cycles
    [max_V,max_ind] = max(avg_voltage);
    [min_V,min_ind] = min(avg_voltage);
    
    % Reorder data by first downwards and then upwards sweep
    grouped_data         = [avg_voltage,I_cycles];
    downwards_sweep_data = grouped_data(max_ind:min_ind-1,:);
    upwards_sweep_data   = [grouped_data(min_ind:end,:); grouped_data(1:max_ind-1,:)];
    half_sweep_ind       = min_ind - max_ind;
    full_sweep_data      = [downwards_sweep_data; upwards_sweep_data];
    
    % Retreive average current, orderded as downwards + upwards sweep
    avg_current          = mean(full_sweep_data(:,2:end),2);    % mean over individual cycles
    std_current          = std(full_sweep_data(:,2:end),[],2);  % std over individual cycles
    avg_voltage          = full_sweep_data(:,1);
    
    % Area calculation
    min_size        = min(size(upwards_sweep_data,1), size(downwards_sweep_data,1));
    areas           = upwards_sweep_data(1:min_size,2:end) - flipud(downwards_sweep_data(1:min_size,2:end));
    flat_loop       = mean(areas, 2);
    flat_loop_std   = std(areas,[],2);
    flat_v          = upwards_sweep_data(1:min_size,1);
    
    % Ratio A+ / A-
    min_voltage_index   = find(min(avg_voltage));
    Aneg_loop           = sum(flat_loop(1:round(end/2))) / numel(flat_loop(1:round(end/2)));
    Apos_loop           = sum(flat_loop(round(end/2)+1:end)) / numel(flat_loop(round(end/2)+1:end));
    ratio               = Apos_loop / Aneg_loop;
    
    figure_title = strrep(figure_title_format,'_',' ');
    figure_sub_title = normalize_by + " sine "+freq+"Hz "+amp+"mV n="+n;
    
    figure()
    hold on
    
    % Plot filled errorbars
    plot_fill(avg_voltage(1:half_sweep_ind),avg_current(1:half_sweep_ind),std_current(1:half_sweep_ind),fill_color)              % positive loop
    plot_fill(avg_voltage(half_sweep_ind+1:end),avg_current(half_sweep_ind+1:end),std_current(half_sweep_ind+1:end),fill_color)  % negative loop
    
    % Plot average
    plot(avg_voltage,avg_current, 'color', 'black','LineWidth',2);
    hold off
        
    set(gca,'fontsize',20)
    xlabel('mV')
    ylabel(y_label) 
    xlim([-amp-2 amp+2])
    title(figure_title);
    subtitle(figure_sub_title);
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    set(gca,'TickDir','both')
    box off
    [~,fig_property] = format_figure(figure_width,figure_height,figure_fontsize,format);
    if export_plots == true
        figure_file_title = strrep(figure_title + " " + figure_sub_title,' ', '_');
        hgexport(gcf,fullfile(plotting_path,figure_file_title+"."+format),fig_property); 
    end
    figure(gcf)
   
    % Plot the individual I cycles
    figure()
    title(figure_title);
    colororder({'k','r'})
    t_0 = 0;
    for w=1:size(I_cycles,2)
        I_cycle = I_cycles(:,w);
        time_to_plot = t_0 + linspace(0,1/freq,length(I_cycle));
        
        yyaxis right
        plot(time_to_plot, av_voltage,'-r','LineWidth',2);
        yyaxis left
        plot(time_to_plot, I_cycles(:,w),'-k','LineWidth',2);
        t_0 = t_0 + 1/freq;
        hold on
    end
    yyaxis left
    ylabel(y_label)
    yyaxis right
    ylabel('mV')
    xlabel('Time (s)')
    
    set(gca,'fontsize',20)
    ax = gca;
    ax.XAxisLocation = 'origin';
    set(gca,'TickDir','both') 
    box off
    
    [~,fig_property] = format_figure(figure_width,figure_height,figure_fontsize,format);
    if export_plots == true
        hgexport(gcf,fullfile(plotting_path,figure_file_title+'  I-traces.'+format),fig_property);
    end
    figure(gcf)
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SQUARE WAVES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Combine up_steps of all channels
up_steps = initialize_dataframe(fieldnames(data(1).up_steps));
for i = 1:length(data)
    up_steps = [up_steps, data(i).up_steps];
end
up_steps = up_steps(2:end); % Remove first empty field

% Combine down_steps of all channels
down_steps = initialize_dataframe(fieldnames(data(1).down_steps));
for i = 1:length(data)
    down_steps = [down_steps, data(i).down_steps];
end
down_steps = down_steps(2:end); % Remove first empty field

%% Do up-steps and down-steps seperately
num_A_points_tested = 1e3; % for exponential fits
step_directions = {'up','down'};
for s = 1:2 % First up, then down
    % 1st iteration deals with up steps, 2nd iteration with down steps
    if s==1
        cycle_data = up_steps;
    else
        cycle_data = down_steps;
    end
    
    % ID vector that assigns each amplitude
    a = [cycle_data(:).amp]';
    [~,~,ID] = unique(a,'rows','stable');

    % Go through every unique amplitude
    for i = 1:length(unique(ID)) % goes through every unique a
        n = 0;                   % sets back number of cycles for next iteration
        unique_a_indeces = find(ID == i); % lists all instances of a given unique combi of a & f
        % Find maximal cycle size (to be used for padding later on):
        max_cycle_size = 0;
        for j = 1:length(unique_a_indeces)
            k = unique_a_indeces(j);
            max_cycle_size = max(max_cycle_size, size(cycle_data(k).V_cycle,1));
        end
        amp               = a(unique_a_indeces(1),1);
        V_cycles          = nan(max_cycle_size, length(unique_a_indeces));   % Store V of each cycle in one column
        I_cycles          = nan(max_cycle_size, length(unique_a_indeces));   % Store I of each cycle in one column
        n_cycles          = nan(max_cycle_size, length(unique_a_indeces));   % Store n of each cycle in one column
        t_cycles          = nan(max_cycle_size, length(unique_a_indeces));   % Store t of each cycle in one column
        cycles            = zeros(max_cycle_size, 3);                        % One column for voltage, one for current, one for time
        
        % If conductance data is available, estimate conductance at this voltage. 
        % Else, use the rough value entered by the user.
        if ~isequal(g_data,'none')
            g_open_est = interp1(g_data(:,1),g_data(:,3),amp*1e-3);
        else
            g_open_est = g_open;
        end
        if isnan(g_open_est)
            disp(['WARNING: conductance data can not be interpolated at V = ',num2str(amp),' mV. Using user-entered open-pore conductance of ',num2str(g_open),' S.'])
            g_open_est = g_open;
        end
        % Loop over cycles with this amplitude
        for j = 1:length(unique_a_indeces)
            k = unique_a_indeces(j);
            cycle = [cycle_data(k).V_cycle, cycle_data(k).I_cycle, cycle_data(k).t_cycle - cycle_data(k).t_cycle(1)];
            sampling_rate = cycle_data(k).sampling_rate;

            % used for padding in case cycle arrays are not same length
            if size(cycle,1) < max_cycle_size
                pad = size(cycles,1) - size(cycle,1);
                cycle = [cycle;nan(pad,3)];
            end
            
            % Estimate number of pores
            num_pores  = 1e-9*max(abs(cycle(:,2))) / (abs(amp)*1e-3) / g_open_est;

            % Do not include this cycle if there were too few pores
            if num_pores < min_num_pores
                disp(['Excluded up-step with num_pores = ',num2str(num_pores), 'S. This number of pores is too low.'])
                continue
            end

            n = n + 1;                              % indicates how many cycles are summed
            cycles = cycles(:,1:3) + cycle(:,1:3);  % sums the values of all cycles
            
            % Compute and save data
            n_cycles(:,n) = (1e-9*cycle(:,2)/(num_pores*amp*1e-3) - g_closed) / (g_open_est - g_closed);
            V_cycles(:,n) = cycle(:,1);
            t_cycles(:,n) = cycle(:,3);
            if isequal(normalize_by,'max') || isequal(normalize_by,'slope')
                I_cycles(:,n) = cycle(:,2) / max(abs(cycle(:,2)));  % normalize the current by max. Normalization by slope does not make sense in the case of a square wave
                y_label = 'norm.I';
            elseif isequal(normalize_by,'none')
                I_cycles(:,n) = cycle(:,2);                         % No normalization
                y_label = 'I (nA)';
            else
                error('Incorrect value for normalize_by. Choose between "max" and "slope"');
            end
            amp = a(k);
        end

        if n==0
            disp("No valid up steps to be plotted for A="+amp)
            continue
        end
        disp("Plotting "+step_directions{s}+"-steps for A="+amp)
        
        % Reduce size of arrays to n, the number of valid cycles
        V_cycles = V_cycles(:,1:n);
        I_cycles = I_cycles(:,1:n);
        n_cycles = n_cycles(:,1:n);
        t_cycles = t_cycles(:,1:n);
        
         % If you want to sum cycles or do no normalization, overwrite I_cycles:
        if isequal(combine_cyles,'sum')  % combine cycles by summing. This weighs cycles with higher current more!
            if isequal(normalize_by, 'max') || isequal(normalize_by,'slope')
                norm_current    = cycles(:,2) / max(abs(cycles(:,2))); % normalizing summed current by its max. Normalizing by slope does not make sense for square waves
            elseif isequal(normalize_by,'none')
                norm_current    = cycles(:,2);
            else
                error('Incorrect value for normalize_by. Choose between "max" and "slope"');
            end
            I_cycles    = norm_current; % sum of the current, normalized
            num_pores   = 1e-9*max(abs(I_cycles)) / (abs(amp)*1e-3) / g_open_est;
            n_cycles    = (1e-9*I_cycles/(num_pores*amp*1e-3) - g_closed) / (g_open_est - g_closed);
        end
        
        % Calculate averages
        av_voltage           = cycles(:,1) / n; % averaging the voltage
%         av_time              = nanmean(t_cycles,2);
%         avg_current          = nanmean(I_cycles,2);    % mean over individual cycles
%         std_current          = nanstd(I_cycles,[],2);  % std over individual cycles
%         avg_n                = nanmean(n_cycles,2);
%         std_n                = nanstd(n_cycles,[],2);  % std over individual cycles
        avg_voltage          = av_voltage;
        
        av_time              = mean(t_cycles,2,'omitnan');
        avg_current          = mean(I_cycles,2,'omitnan');    % mean over individual cycles
        std_current          = std(I_cycles,[],2,'omitnan');  % std over individual cycles
        avg_n                = mean(n_cycles,2,'omitnan');
        std_n                = std(n_cycles,[],2,'omitnan');  % std over individual cycles

        % Plot filled errorbars
        figure_title = figure_title_format + " "+step_directions{s}+"-step "+amp+"mV n="+n;
        figure()
        hold on
        
        % Fit exponential:
        % If conductance data is available, then plot the estimated number
        % of pores and the exponential fit on it.
        % If not, then plot the normalized current and the exponential fit
        % on it.
        if ~isequal(g_data,'none')
            % Plot estimated nr of pores
            notnan  = ~isnan(avg_n);
            plot_fill(av_time(notnan),avg_n(notnan),std_n(notnan),fill_color)
            p1 = plot(av_time(notnan),avg_n(notnan), 'color', 'black','LineWidth',2);
            % Fit exponential decay on estimated number of pores
            [A,tau] = estimate_tau_on_log_scale(avg_n(notnan),av_time(notnan),1,num_A_points_tested);
            r_open  = A/tau;
            r_close = (1-A) / tau;
            leg = {'Data',['2-state model fit (r_o=',num2str(r_open,2),', r_c=',num2str(r_close,2),')']};
            ylabel('Estimated fraction of open pores')
        else
            % Plot avg current
            notnan  = ~isnan(avg_current);
            plot_fill(av_time(notnan),avg_current(notnan),std_current(notnan),fill_color)
            p1 = plot(av_time(notnan),avg_current(notnan), 'color', 'black','LineWidth',2);
            % Fit exponential decay on current normalized by max
            notnan  = ~isnan(avg_current);
            [A,tau] = estimate_tau_on_log_scale(avg_current(notnan),av_time(notnan),1,num_A_points_tested);
            A = A - 2*(s==2)*A; % flip the sign of A if this is a downwards step
            leg = {'Data',['exponential fit (tau=',num2str(tau,2),', A=',num2str(A,2),')']};
            ylabel(y_label)
        end
        t = linspace(av_time(1),av_time(end),100);
        if s==1 % up-step
            p2 = plot(t,A+(1-A)*exp(-t/tau),'color','cyan','LineWidth',2);
        else % down-step
            p2 = plot(t,A-(1+A)*exp(-t/tau),'color','cyan','LineWidth',2);
        end
        legend([p1,p2],leg,'Location','best')
        hold off

        % Figure formatting and export
        ylabel('Estimated fraction of open pores')
        xlabel('Time (s)')
        title(strrep(figure_title,'_', ' '));
        set(gca,'TickDir','both')
        box off
        [~,fig_property] = format_figure(figure_width,figure_height,figure_fontsize,format);
        if export_plots == true
            figure_file_title = strrep(figure_title,' ', '_');
            hgexport(gcf,fullfile(plotting_path,figure_file_title+"."+format),fig_property); 
        end
        figure(gcf)

        % Plot the individual I cycles
        figure()
        title(strrep(figure_title,'_', ' '));
        colororder({'k','r'})
        t_0 = 0;
        for w=1:size(I_cycles,2)
            I_cycle = I_cycles(:,w);
            t_frame = (0:length(I_cycle)-1)/sampling_rate;
            time_to_plot = t_0 + t_frame;
            yyaxis right
            plot(time_to_plot, av_voltage,'-r','LineWidth',2);
            yyaxis left
            plot(time_to_plot, I_cycles(:,w),'-k','LineWidth',2);
            t_0 = t_0 + t_frame(end);
            hold on
        end
        
        % Figure formatting and export
        yyaxis left
        ylabel(y_label)
        yyaxis right
        ylabel('mV')
        xlabel('Time (s)')
        set(gca,'TickDir','both') 
        box off
        [~,fig_property] = format_figure(figure_width,figure_height,figure_fontsize,format);
        if export_plots == true
            hgexport(gcf,fullfile(plotting_path,figure_file_title+'  I-traces.'+format),fig_property);
        end
        figure(gcf)
    end
end