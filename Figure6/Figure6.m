%% Functions to plot figure 6
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
    parameters.total_num_animals = 13;
    parameters.new_region_idx = [1 2 3 5 6 4 7 8];
end

%% Plot enrichment index
load 'matched_dataset.mat'

% Get animal-averaged enrichment index across all regions
enrichment = get_enrichment_idx(parameters,data,data.cell_idx.coupling_sig_cells);

% Plot scatterplot for enrichment index
plot_scatter_enrich_idx(enrichment.averaged_sensory) % Activity projected on stimulus mode
plot_scatter_enrich_idx(enrichment.averaged_choice) % Activity projected on choice mode

% Plot AUC correlation
plot_AUC_correlation(data,enrichment)

% Get regional enrichment index
regional_enrichment = get_regional_enrichment_idx(parameters,data,data.cell_idx.coupling_sig_cells);

% Plot regional enrichment index heatmap
plot_regional_enrich_idx(parameters,regional_enrichment.stim_right) % For right-stimulus encoding cells
plot_regional_enrich_idx(parameters,regional_enrichment.delay_right) % For right-delay encoding cells

% Plot regional enrichment index network
plot_regional_idx_network(parameters,regional_enrichment.stim_right) % For right-stimulus encoding cells
plot_regional_idx_network(parameters,regional_enrichment.delay_right) % For right-delay encoding cells

%% Plot ablated population activity
% Get reconstructed ablated data and project onto choice axes (Long runtime)
cd([base,'/Ablation/Raw']);
load 'ablated_data_raw.mat'
ablated_data = get_ablated_data(parameters,data);

% Or Load processed dataset
cd([base,'/Ablation'])
load 'ablated_data.mat'

% Plot choice mode after selective ablation
% Selected regions ALM, vS1, RSC and PPC
plot_ablated_choice_mode(parameters,ablated_data)

% Plot choice epoch-averaged choice mode
% Bar and heatmap
plot_ablated_epoch_averaged_choice_mode(parameters,ablated_data)