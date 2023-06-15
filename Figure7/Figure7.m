%% Functions to plot figure 7
% Initialize parameters:
for hide = 1
    parameters.root_folder = uigetdir; % CD to folder
    cd(parameters.root_folder);

    % Experiment parameters
    parameters.fs_image = 9.3521;
    parameters.pre_resp_epoch = [1:round(parameters.fs_image*4)];
    parameters.post_sample_epoch = [round(parameters.fs_image*2):round(parameters.fs_image*4)];
    parameters.choice_epoch = [round(parameters.fs_image*3):round(parameters.fs_image*4)];
    parameters.sample_epoch = round(parameters.fs_image*1) : round(parameters.fs_image*2);
    parameters.ITI_epoch = 1:round(parameters.fs_image);
    parameters.sample_offset = round(parameters.fs_image*2);
    parameters.action_epoch = [round(parameters.fs_image*4):round(parameters.fs_image*5)];
    parameters.baseline_epoch = parameters.ITI_epoch;
    parameters.norm_epoch = parameters.pre_resp_epoch;
    parameters.total_num_animals = 13;
    parameters.new_region_idx = [1 2 3 5 6 4 7 8];
end

%% Optogenetics
load optogenetics_data.mat
[ppc_correct_rate, vs1_correct_rate] = analyze_behaviour_opto(opto_data);

% Plot trained behaviour data 
session_idx = [1 2 3]; % First 3 sessions without overtraining
opto_behaviour = plot_behaviour_opto(ppc_correct_rate,vs1_correct_rate,session_idx);

% Plot overtrained behaviour data 
session_idx = [4 5 6]; % Last 3 sessions after overtraining
opto_behaviour = plot_behaviour_opto(ppc_correct_rate,vs1_correct_rate,session_idx);

% Plot bar plot for effects of overtraining for PPC
plot_diff_overtrained(ppc_correct_rate)

%% DREADD
load 'DREADD_data.mat'
% Get DREADD behaviour
[CNO_correct_rate, Saline_correct_rate] = analyze_behaviour_dreadd(DREADD_data);

% Plot DREADD performance
[concat_data,p_value] = plot_behaviour_DREADD(CNO_correct_rate,Saline_correct_rate);

% Get task-encoding neurons
load 'DREADD_activity.mat'
% Get sample, choice and action encoding cells
task_encoding_cells = get_DREADD_task_encoding_neurons(parameters,DREADD_activity);

% Get activity for choice-encoding cells
average_activity = get_mean_activity_choice_encoding(parameters,DREADD_activity,task_encoding_cells);