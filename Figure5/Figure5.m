%% Functions to plot figure 5
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
    parameters.new_region_idx = [1 2 3 5 6 4 7 8];
end

%% Reconstruct retained or eliminated activity
% Load raw dataset containing raw activity
load 'retained_eliminated_data_raw.mat' 

% Get reconstructed retained and eliminated activity (long runtime)
retained_elim_data = get_retained_eliminated_activity(parameters,data);

% Or load pre-processed dataset
load 'retained_elim_data.mat'

% Plot average number of retained and lost functional coupling
[prop_retained,prop_lost] = plot_num_retained_lost_cells(retained_elim_data);

% Get Choice mode from reconstructed activity
retained_elim_data = get_choice_mode(parameters,data,retained_elim_data);

% Plot averaged choice mode
plot_choice_mode(parameters,retained_elim_data)

%% Get AUC for each shuffled data (long runtime)
decoding_accuracy = get_decoding_accuracy(parameters,retained_elim_data);

% Or load pre-processed data
load('decoding_accuracy.mat')

% Plot difference in decoding accuracy (retained - eliminated)
plot_decoding_accuracy(parameters,decoding_accuracy)