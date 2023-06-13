%% Functions to plot figure 2
% Initialize parameters:
for hide = 1
    base = uigetdir; % CD to folder
    parameters.root_folder =  [base,'/Non_matched'];
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

%% Extract data from raw dataset and project onto sensory and choice axes (long runtime)
data = get_data(parameters);

% Or directly load the processed dataset (with choice and sensory mode)
load('non_matched_dataset.mat'); % Choice -> stim orthogonalisation (default)

%% Get ROC for each animal (Population level)
% Get AUC using ROC analysis for stimulus/choice population activity projected on stimulus or choice axes.
% 'choice' Decode based on animal's choice 
% 'stimulus' Decode based on trial type
ROC_population_data.stimulus = get_ROC_choice_stim_mode(data,'stimulus');
ROC_population_data.choice = get_ROC_choice_stim_mode(data,'choice');

% Plot decoded accuracy
plot_auc_choice_stim_population(parameters,ROC_population_data.stimulus,parameters.sample_epoch)
plot_auc_choice_stim_population(parameters,ROC_population_data.choice,parameters.post_sample_epoch)