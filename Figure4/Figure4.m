%% Functions to plot figure 4
% Initialize parameters:
for hide = 1
    base = uigetdir; % CD to folder
    parameters.root_folder =  [base,'/Matched'];
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

%% Analyse GLM weights
load 'matched_dataset.mat'

% Get GLM weights to analyse
% 'stimulus/delay/action' - Task variable predictors
GLM_coeff = get_GLM_coeff(parameters,data,'stimulus');
GLM_coeff = get_GLM_coeff(parameters,data,'delay');
GLM_coeff = get_GLM_coeff(parameters,data,'action');

%% Analyse convergence index
% Get conv idx data
convergence_idx = get_conv_idx(data,data.cell_idx.coupling_sig_cells);

% Plot graphs for global conv. idx
plot_global_conv_idx(convergence_idx)

% Plot graphs for regional conv. idx
plot_conv_idx_heatmap(parameters,convergence_idx)

% Plot graphs for conv. idx network
plot_conv_idx_network(parameters,convergence_idx)