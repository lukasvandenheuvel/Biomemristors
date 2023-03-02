clear all
close all
addpath('./functions')

%% Intro
% Hello! This function will search for (1) sine waves, (2) upward steps and
% (3) downward steps in the voltage traces you recorded.
% You will go through them and accept or trash them.
% The accepted cycles will be saved, together with the original current and
% voltage traces.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% SET PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

input_path      = '/Volumes/lben-archives/2023/Simon/OrbitMiniData';
output_path     = '/Volumes/lben-archives/2023/Simon/Biomemristors/Biomemdata-checked/AelK238AK242A';
%filename        = 'AelK238AK242A_KCl_1M_pH6p2_10C_100us__230202131631';
filename        = 'AelK238AK242A_KCl_1M_pH6p2_20C_100us__230202125956';
filename_format = 'PORE_SALT_CONCENTRATION_PH_TEMPERATURE_PERIOD_EMPTY_TIMESTAMP'; % how the filename is built up (seperate snippets with '_'). For femto data, a field PERIOD must exist which specifies the sampling period.
device          = 'orbitmini';              % choose 'femto', 'axopatch' or 'orbitmini'
figure_position = [0,0,0.9,0.9];            % figure position: [x,y,width,height]. x and y are positions of the lower-left figure corner. Values must be between 0 and 1, relative to the screen size. 

clipping_length_threshold   = 10;   % regions where the current doesn't change for longer than clipping_length_threshold indeces are considered clipping events.
fft_peak_theshold           = 0.75; % Decrease this if not all sine waves are properly found.
min_step_size               = 9;    % minimal size of voltage step in mV
min_step_duration           = 0.5;  % minimal duration of a voltage step in sec
membrane_relax_time         = 50;   % after a voltage step, the time it takes before the membrane is charged in ms.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% START SCRIPT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Init -- check input and output folders
% Check if the input folder exists
assert(isfolder(input_path), 'The input path could not be found!');
% Check if the output folder exists. If not, create it.
if ~isfolder(output_path)
    mkdir(output_path);
    disp(['Created a new folder ',output_path])
end
[~, filename_no_ext, ~] = fileparts(filename);
output_filepath         = fullfile(output_path, [filename_no_ext,'.mat']);

%% Read data
disp("Loading data ...")
input_filepath  = fullfile(input_path, filename);
data            = read_data(input_filepath,device,filename_format);
if isempty(data(1).V)
    disp(['No voltage or current trace was found in ',input_filepath,'. Aborting script.'])
    return
end

%% Init
data_clean_fields   = {'V','I','sampling_rate','header','conditions','device','sine_waves','up_steps','down_steps'};
data_clean          = initialize_dataframe(data_clean_fields);
cycle_types_lut     = {'sine','up','down'};                                 % Three cycle types are possible: sine wave, up-step, down-step
cycle_data_fields   = {'freq','amp','duration','type','I_cycle','V_cycle','t_cycle','sampling_rate'}; % fields to save for each cycle

%% Obtain screen resolution
set(0,'units','pixels')         % Sets the units of your root object (screen) to pixels
pix_ss = get(0,'screensize');   % Obtains this pixel information

assert(~any(figure_position<0) & ~any(figure_position>1), 'figure_position must contain values between 0 and 1, relative to screen size.')
figure_x      = round(pix_ss(3) * figure_position(1));
figure_y      = round(pix_ss(4) * figure_position(2));
figure_width  = round(pix_ss(3) * figure_position(3));
figure_height = round(pix_ss(4) * figure_position(4));

%% Display instructions
disp('....................................................................................')
disp('Hello, with this program you will be able to select cycles to include in later analysis.')
disp(['The chosen waves will be saved in ',output_filepath])
disp('You will select cycles with keys on your keyboard:')
disp('.................................')
disp('Right arrow → = include this cycle')
disp('Down arrow  ↓ = exclude this cycle')
disp('Left arrow  ← = go back to previous')
disp('Enter         = Leave program. The data will be saved.')
disp('.................................')
disp('To show this message again, press any other random key.')

%% Loop through channels
for c = 1:length(data)
    disp(['Opening trace of channel ',num2str(c),'. Please wait ...'])
    % data_clean inherits all properties of data
    data_clean(c).V             = data(c).V;
    data_clean(c).I             = data(c).I;
    data_clean(c).sampling_rate = data(c).sampling_rate;
    data_clean(c).header        = data(c).header;
    data_clean(c).conditions    = data(c).conditions;
    data_clean(c).device        = data(c).device;

    % Init variables which depend on sampling rate and exp. conditions
    sampling_period = 1/data(c).sampling_rate;
    n_relax         = round(1e-3*membrane_relax_time * data(c).sampling_rate);
    t               = sampling_period*(0:length(data(c).I)-1)';
    W_voltage_mean  = round(5 / sampling_period / 100);      % window size for averaging
    fig_title       = strjoin(struct2cell(data(c).conditions));

    % Remove current clipping regions (i.e., regions where membrane breaks) 
    clipped = find_clipped_regions(data(c).I, clipping_length_threshold);
    I = data(c).I(~clipped);
    V = data(c).V(~clipped);
    t = t(~clipped);

    % Find regions of valid sine waves
    sine_wave_inds = find_sine_waves(movmean(V,W_voltage_mean),fft_peak_theshold); % An Nx2 matrix with the start and end index of each identified wave

    % Find regions of square waves (up and down steps)
    min_step_length = round(min_step_duration * data(c).sampling_rate);         % minimal duration of a voltage step (indeces)
    up_step_inds    = find_up_steps(V, min_step_size, min_step_length);         
    down_step_inds  = find_down_steps(V, min_step_size, min_step_length);

    % Combine all waves into one array
    cycle_inds_unsorted  = [sine_wave_inds; up_step_inds; down_step_inds];
    cycle_types_unsorted = [zeros(size(sine_wave_inds,1),1) + 1; ...
                            zeros(size(up_step_inds,1),1) + 2; ...
                            zeros(size(down_step_inds,1),1) + 3];
                        
    if isempty(cycle_inds_unsorted)
        disp(['No cycles found in ',fig_title,' channel ',num2str(c)])
        continue
    end

    [~,inds]    = sort(cycle_inds_unsorted(:,1));
    cycle_inds  = cycle_inds_unsorted(inds,:);
    cycle_types = cycle_types_unsorted(inds,:);
    num_waves   = size(cycle_inds,1);
    
    % Initialize empty data structure
    empty_data  = cell(length(cycle_data_fields),1);
    cycle_data  = cell2struct(empty_data, cycle_data_fields);
    
    % Initialize other variables
    answer_history = zeros(1,num_waves);
    cycles_in_data = [];
    count = 0;  % keeps track of nr of cycles still to go
    n = 0;      % keeps track of nr of cycles in output data

    % Init plot
    figure('color','w','units','points','position',[figure_x,figure_y,figure_width,figure_height])

    subplot(2,3,[1,2])
    plot(t,I,'k')
    hold on
    p1 = plot(t(1),I(1),'r','LineWidth',2);
    ylabel('Current (nA)')
    xlabel('Time (s)')

    subplot(2,3,[4,5])
    plot(t,V,'k')
    hold on
    p2 = plot(t(1),V(1),'r','LineWidth',2);
    ylabel('Voltage (mV)')
    xlabel('Time (s)')

    subplot(2,3,3);
    p3 = plot(0,0, 'k');
    ylabel('Current (nA)')
    xlabel('Time (s)')

    subplot(2,3,6);
    p4 = plot(0,0,'k');
    ylabel('Current (nA)')
    xlabel('Voltage (mV)')

    skipped_channel = false;

    while count < num_waves % continue until all waves are done
        % Increment
        j = j + 1;
        count = count + 1;
        % Obtain V, I and t in cycle
        cycle_type  = cycle_types(count);
        is_square   = ((cycle_type == 2) || (cycle_type == 3));
        start_index = cycle_inds(count,1) + is_square*n_relax;              % Shorten square waves by the membrane relaxation time
        end_index   = cycle_inds(count,2) - is_square*n_relax;              % Shorten square waves by the membrane relaxation time
        V_cycle = V(start_index:end_index);
        I_cycle = I(start_index:end_index);
        t_cycle = t(start_index:end_index);
        % Obtain frequency and amplitude
        freq    = round(1/((end_index - start_index)*sampling_period),1,'significant');     % frequency
        duration = 1/freq;
        if cycle_type==3 % down step
            amp = 100*round(min(V_cycle)/100,1); % amplitude, rounded by tens
        else
            amp = 100*round(max(V_cycle)/100,1); % amplitude, rounded by tens
        end
        time    = linspace(0,1/freq,length(V_cycle));

        % Update plots
        p1.XData = t_cycle;
        p1.YData = I_cycle;
        p2.XData = t_cycle;
        p2.YData = V_cycle;
        p3.XData = time; 
        p3.YData = I_cycle;
        p4.XData = V_cycle;
        p4.YData = I_cycle;
        suptitle({fig_title, ['Channel ',num2str(c),' / ',num2str(length(data))],['Cycle ',num2str(count),' / ',num2str(num_waves),' = ',cycle_types_lut{cycle_type}]})

        % Wait for keyboard press
        press  = waitforbuttonpress;
        answer = double(get(gcf,'CurrentCharacter'));
        answer_history(count) = answer;

        switch answer
            
            case 115 % skip channel = 's'
                disp(['Channel ',num2str(c),' was skipped.'])
                skipped_channel = true;
            
            case 29 % add = right arrow
                n = n + 1;  % indicates how many cycles are added to the structure
                cycle_data(n).sampling_rate = data(c).sampling_rate;
                cycle_data(n).freq      = freq;
                cycle_data(n).amp       = amp;
                cycle_data(n).duration  = duration;
                cycle_data(n).type      = cycle_type;
                cycle_data(n).I_cycle   = I_cycle;
                cycle_data(n).V_cycle   = V_cycle;
                cycle_data(n).t_cycle   = t_cycle;
                cycles_in_data(n)       = count;
                % save(output_filepath,'cycle_data')
                disp(['Cycle ',num2str(count),' was added.'])
            case 31 % trash = down arrow
                disp(['Cycle ',num2str(count),' was trashed.'])
            case 28 % go back = left arrow
                count = count - 2;
                % if the cycle was added in the previous round, then remove it now
                if answer_history(count+1)==29
                    n = n - 1;
                    cycle_data     = cycle_data(1:n);
                    cycles_in_data = cycles_in_data(1:n);
                end
                answer_history(count+1) = 0;
                disp(['Go back to cycle ',num2str(count+1)])
            case 13 % cancel program = enter
                disp('Saving data ...')
                data_clean(c).sine_waves = cycle_data([cycle_data(:).type]==1);
                data_clean(c).up_steps   = cycle_data([cycle_data(:).type]==2);
                data_clean(c).down_steps = cycle_data([cycle_data(:).type]==3);
                save(output_filepath,'data_clean')
                disp(['You cancelled the program. The selected waves were saved in ',output_filepath])
                close all;
                return
            otherwise
                j = j-1;
                count = count-1;
                disp('You pressed an invalid key! These are the valid ones:')
                disp('.................................')
                disp('Right arrow → = include this wave')
                disp('Down arrow  ↓ = exclude this wave')
                disp('Left arrow  ← = go back to previous')
                disp('Enter         = Leave program. The data will be saved.')
                disp('.................................')
                disp('To show this message again, press any other random key.')
        end
        
        if skipped_channel
            break
        end
    end
    % Save cycles
    data_clean(c).sine_waves = cycle_data([cycle_data(:).type]==1);
    data_clean(c).up_steps   = cycle_data([cycle_data(:).type]==2);
    data_clean(c).down_steps = cycle_data([cycle_data(:).type]==3);
    close();
end

disp('Saving data ... Please be patient.')
save(output_filepath,'data_clean')
disp(['Done! All selected waves are saved in ',output_filepath])
