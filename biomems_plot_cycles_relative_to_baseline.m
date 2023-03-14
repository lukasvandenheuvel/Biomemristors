clear all
close all

%% Intro
% Hello! With this function you will plot IV-curves normalized by their
% open-pore baseline.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% SET PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

input_path              = '/Users/lukasvandenheuvel/Documents/EPFL/MA3/LBM/Biomemdata_cleaned/Ael-K238AK242A';
iv_path                 = '/Users/lukasvandenheuvel/Documents/EPFL/MA3/LBM/Biomemdata_cleaned/gating&IV/AelwtdoubleKtoA/AelwtdoubleKtoA.csv';  % Enter a path to a csv file
plotting_path           = '/Users/lukasvandenheuvel/Documents/EPFL/MA3/LBM/Biomemdata_plots';
filename_includes       = 'sine_f0.1Hz_C0_M1_Cm52.8_16-01-2023_17-01-06';          % All filenames with this string included will be combined
output_title            = 'filename_format';     % specify a title for output plots OR enter 'filename_format' to use the format of the first loaded datafile.

% analysis parameters
num_cycles              = 15;       % number of cycles to include in analysis. The cycles with the highest maximal current will be picked.
min_num_pores           = 1;        % minimal (estimated) number of pores needed to take measurement into account
V_window_no_gating      = [-35,35]; % Voltage window (in mV) in which almost no gating takes place. Is used to estimate number of pores by fitting the IV-curve in this window to the IV curve of the data.

% figure parameters
format                  =  'pdf';                   % figure extension. png takes a long time
export_plots            =  true;                    % if trdaue plots are saved in the specified format
figure_width            = 15; % inches
figure_height           = 10; % inches
figure_fontsize         = 11;
fill_color              = [0.9, 0.9, 0.9]; % color for errorbars

%% Load data
data = read_mat(input_path, filename_includes);
if isempty(data)
    disp(['No data found folder ',input_path,' corresponding to the name ',filename_includes,'!'])
    return
end

%% Load IV-curve
iv_curve = readtable(iv_path);

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
figure_title_format = figure_title_format;
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
        
        % Estimate I_offset and number of pores by fitting the current to
        % the IV-curve in a voltage window where no gating takes place.
        [I_offset, num_pores_estimate] = fit_current_to_baseline(cycle(:,2), cycle(:,1), iv_curve, V_window_no_gating);
        cycle(:,2) = cycle(:,2) - I_offset;                                 % Remove offset from current
        
        % used for padding in case cycle arrays are not same length
        if size(cycle,1) < max_cycle_size
            pad = size(cycles,1) - size(cycle,1);
            cycle = [cycle;zeros(pad,2)];
        end
        
        % Do not include this cycle if there were too few pores
        if num_pores_estimate < min_num_pores
            disp(['Excluded cycle with estimated N = ',num2str(num_pores_estimate), ' pores. This is not enough pores.'])
            continue
        end
               
        % Save data
        n = n + 1;                                                          % indicates how many cycles are summed
        cycles = cycles(:,1:2) + cycle(:,1:2);                              % sums the values of all cycles
        V_cycles(:,n) = cycle(:,1);
        I_cycles(:,n) = cycle(:,2) / num_pores_estimate;                    % Normalize current by the estimated number of pores 
    end
    
    V_cycles = V_cycles(:,1:n);
    I_cycles = I_cycles(:,1:n);
    
    if n==0
        disp("No valid waves to be plotted for A="+amp+" and f="+freq)
        continue
    end
    disp("Plotting sine waves for A="+amp+" and f="+freq)
    
    % Average data over cycles
    avg_voltage     = mean(V_cycles,2);   % mean over individual cycles
    avg_current     = mean(I_cycles,2);   % mean over individual cycles
    [max_V,max_ind] = max(avg_voltage);
    [min_V,min_ind] = min(avg_voltage);
    
    % Interpolate the IV-curve to make a baseline
    I_baseline = interp1(iv_curve.V_axis_mV, iv_curve.I_pore_nA, avg_voltage);
    
    % Divide cycle in 4 parts
    grouped_data            = [avg_voltage,avg_current, I_baseline];
    [I1,I2,I3,I4]           = split_cycle(avg_current, avg_voltage);
    [Ibl1,Ibl2,Ibl3,Ibl4]   = split_cycle(I_baseline, avg_voltage);
    [V1,V2,V3,V4]           = split_cycle(avg_voltage, avg_voltage);

    Ibc1 = Ibl1 - I1; % baseline - current from V=0 to V=max
    Ibc2 = Ibl2 - I2; % baseline - current from V=max to V=0
    Ibc3 = Ibl3 - I3; % baseline - current from V=0 to V=min
    Ibc4 = Ibl4 - I4; % baseline - current from V=min to V=0
    
    % Plot the baseline-corrected current traces
    figure_title = strrep(figure_title_format,'_',' ');
    figure_sub_title = "sine "+freq+"Hz "+amp+"mV n="+n;
    
    % IV-curve
    figure()
    p1 = plot([V3;V1], [I3;I1], '.k', 'LineWidth', 2);
    hold on
    p2 = plot([V4;V2], [I4;I2], '.b', 'LineWidth', 2);
    p3 = plot(iv_curve.V_axis_mV, iv_curve.I_pore_nA, 'r', 'LineWidth', 2);
    ylabel('norm.I (nA)')
    xlabel('mV')
    title([figure_title, figure_sub_title, ' IV-curve'])
    legend([p1,p3],['Average of ',num2str(n),' pores'], 'Open pore', 'Location', 'NorthWest')
    [~,fig_property] = format_figure(figure_width,figure_height,figure_fontsize,format);
    if export_plots == true
        figure_file_title = strrep(figure_title,' ', '_');
        hgexport(gcf,fullfile(plotting_path,figure_file_title+" VI."+format),fig_property); 
    end
    figure(gcf)

    % Closing 
    figure()
    plot([V3;V1],[Ibc3;Ibc1], '.k')
    title([figure_title, figure_sub_title, ' Closing'])
    ylabel('norm.I - baseline')
    xlabel('mV')
    [~,fig_property] = format_figure(figure_width,figure_height,figure_fontsize,format);
    if export_plots == true
        figure_file_title = strrep(figure_title,' ', '_');
        hgexport(gcf,fullfile(plotting_path,figure_file_title+" closing."+format),fig_property); 
    end
    figure(gcf)

    % Opening 
    figure()
    plot([V4;V2],[Ibc4;Ibc2], '.b')
    title([figure_title, figure_sub_title, ' Opening'])
    ylabel('norm.I - baseline')
    xlabel('mV')
    [~,fig_property] = format_figure(figure_width,figure_height,figure_fontsize,format);
    if export_plots == true
        figure_file_title = strrep(figure_title,' ', '_');
        hgexport(gcf,fullfile(plotting_path,figure_file_title+" opening."+format),fig_property); 
    end
    figure(gcf)
end
