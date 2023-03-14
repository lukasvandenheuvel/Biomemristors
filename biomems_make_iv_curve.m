clear all
close all

input_path              = '/Users/lukasvandenheuvel/Documents/EPFL/MA3/LBM/Data-checked/gating&IV/AelwtdoubleKtoA';
filename_includes       = 'choose_name_07-03-2023_13-43-55';        % All filenames with this string included will be combined
output_title            = 'AelwtdoubleKtoA.csv';                    % specify a title for output plots OR enter 'filename_format' to use the format of the first loaded datafile.

% IMPORTANT IV-curve making parameters
num_pores               = 4;        % The EXACT number of pores throughout the whole recording. The number of pores may NOT change! Very important to get this right!!!!!
N_window                = 20;       % Window size for median filtering.
add_to_amp              = 1;        % voltage (in mV) to add to the amplitude. This should be 0 unless you have a very good reason to add something (e.g., you forgot to record your IV-curves from -205 to 205 mV instead of from -200 mV to 200 mV).
num_V_points            = 403;      % The number of voltage points on which the V-I curve will be defined. Make sure this number is uneven so that V=0 is a point on the line.

%% Load data
data = read_mat(input_path, filename_includes);
if isempty(data)
    disp(['No data found folder ',input_path,' corresponding to the name ',filename_includes,'!'])
    return
end

%% Combine sine_waves of all channels
sine_waves = initialize_dataframe(fieldnames(data(1).sine_waves));
for i = 1:length(data)
    sine_waves = [sine_waves, data(i).sine_waves];
end
sine_waves = sine_waves(2:end); % Remove first empty field

%% Create I-V curve
amp     = unique([sine_waves(:).amp])+add_to_amp;
freq    = unique([sine_waves(:).freq]);
assert((length(amp)==1 && length(freq)==1), 'The amplitude and/or frequency of the waves changes throughout the recording! Run biomems_check_cycles again to avoid this.');

% Initialize V axis and average I arrays
V_axis_mV   = linspace(-amp,amp,num_V_points)';
I_avg       = zeros(num_V_points, length(sine_waves));

for i = 1:length(sine_waves)
    I = movmedian( sine_waves(i).I_cycle, N_window ) / num_pores;           % current normalized by the number of pores
    V = sine_waves(i).V_cycle;                                              % voltage
    
    % Divide cycle in 4 parts
    [max_V,max_ind] = max(V);
    [min_V,min_ind] = min(V);
    half_ind = round((min_ind+max_ind)/2);
    I1 = I(1:max_ind);                                                      % From V=0 to V=max
    I2 = I(max_ind+1:half_ind);                                             % From V=max to V=0
    I3 = I(half_ind+1:min_ind);                                             % From V=0 to V=min
    I4 = I(min_ind+1:end);                                                  % From V=min to V=0
    
    % Rearrange current into up-cycle and down-cycle
    I_up    = [I4; I1];
    I_down  = [I2; I3];
    V_up    = amp*linspace(-1,1,length(I_up));                              % The voltage is made artificially to make sure there are no unique values. This is important for the interpolation function.
    V_down  = amp*linspace(1,-1,length(I_down));                            % The voltage is made artificially to make sure there are no unique values. This is important for the interpolation function.
    
    % Interpolate voltage
    I_up_interp     = interp1(V_up,I_up,V_axis_mV);                            % Interpolate the current in the up-cycle on the points on V_axis
    I_down_interp   = interp1(V_down,I_down,V_axis_mV);                        % Interpolate the current in the down-cycle on the points on V_axis
    I_avg(:,i)      = (I_up_interp + I_down_interp) / 2;                    % Average up and down cycles to get rid of capacitance
end

% Mean current 
I_total = mean(I_avg,2); % Take mean over cols

% Remove offset
offset      = I_total(V_axis_mV==0);
I_pore_nA   = I_total - offset;

%% Plot I-V curve
figure()
plot_fill(V_axis_mV,I_pore_nA,std(I_avg,[],2),[0.5,0.5,0.5]);
hold on
plot(V_axis_mV,I_pore_nA,'k')
title(['Average I-V curve for ',num2str(num_pores),' pores'])
ylabel('I (nA)')
xlabel('V (mV)')

%% Save I-V curve
T = table(V_axis_mV, I_pore_nA);
out_dir = fullfile(input_path, output_title);
writetable(T,out_dir,'Delimiter',',');
disp(['I-V curve saved to ', out_dir])