%% Functions to plot figure 1
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
    parameters.naive_session = [1 2];
    parameters.inter_session = [3 4];
    parameters.expert_session = [5 6];
    parameters.new_region_idx = [1 2 3 5 6 4 7 8];
end

%% Get ROC for each animal (Single cell level)
load 'matched_dataset.mat'

% Get AUC using ROC analysis for single cell activity.
% 'choice' Decode based on animal's choice 
% 'stimulus' Decode based on trial type
ROC_data.stimulus = get_ROC_stimulus_activity(data,'stimulus');
ROC_data.choice = get_ROC_stimulus_activity(data,'choice');

plot_auc_data(parameters,ROC_data.stimulus,parameters.sample_epoch)
plot_auc_data(parameters,ROC_data.choice,parameters.post_sample_epoch)