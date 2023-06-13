%% Functions to plot figure 3
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

%% Get coordination between pairs of regions
load('non_matched_dataset.mat');

% For both choice and sensory data
coordination_data = get_coordination_data(parameters,data);

%% Plot coordination heatmap
% 'sensory/choice' - coordination between stimulus trials or animal's choice, respectively
% 'same/both' - Within same trial types or across both trial types, respectively
% 'naive/intermediate/expert' - learning stage
% E.g. plot_coordination_heatmap(parameters, data, 'sensory/choice', 'both/same','naive/intermediate/expert','naive/intermediate/expert')
plot_coordination_heatmap(parameters,coordination_data,'sensory','same','naive','expert')
plot_coordination_heatmap(parameters,coordination_data,'choice','same','naive','expert')

plot_coordination_heatmap(parameters,coordination_data,'sensory','both','naive','expert')
plot_coordination_heatmap(parameters,coordination_data,'choice','both','naive','expert')

% Plot p-value network
plot_coordination_network_p_value(parameters,coordination_data,'sensory','same','naive','expert')
plot_coordination_network_p_value(parameters,coordination_data,'sensory','both','naive','expert')

plot_coordination_network_p_value(parameters,coordination_data,'choice','same','naive','expert')
plot_coordination_network_p_value(parameters,coordination_data,'choice','both','naive','expert')