%% (Demo) Coupling analysis combine all animals to get the coupling parameters
% Load demo data
parameters.root_folder = uigetdir;
cd([parameters.root_folder,'/GLM_demo'])
load GLM_parameters.mat
load GLM_result.mat

fs_image = 9.352;
p_value_threshold = 0.001;
no_partial_counter = 1;
session_num = 3; % For demo data
animal_num = 7; % For demo data

% Check if the number of cells is correct for all files
assert(size(GLM_result.coeff,1) == size(resp_task1_var_t1,2),'number of cells is wrong')
assert(size(p_value,2) == size(resp_task1_var_t1,2),'number of cells is wrong')
assert(size(p_value_task1_var,1) == size(resp_task1_var_t1,2),'number of cells is wrong')
assert(size(explained_deviance_all_reduction,1) ==  size(resp_task1_var_t1,2),'number of cells is wrong')

cell_coupling = predictor_labels(~isnan(predictor_labels));
coeff{animal_num}{session_num} = GLM_result.coeff;
all_coeff{animal_num}{session_num} = GLM_result.coeff;

total_cell_num = size(GLM_result.coeff,1);

tact_stim_onset_left_idx = 1:stim_nBases;
tact_stim_onset_right_idx = (max(tact_stim_onset_left_idx)+1) : (max(tact_stim_onset_left_idx)+stim_nBases);
delay_left_middle_idx = (max(tact_stim_onset_right_idx)+1) : (max(tact_stim_onset_right_idx)+delay_nBases);
delay_right_middle_idx = (max(delay_left_middle_idx)+1) : (max(delay_left_middle_idx)+delay_nBases);
left_lick_onset_idx = (max(delay_right_middle_idx)+1) : (max(delay_right_middle_idx)+resp_nBases);
right_lick_onset_idx = (max(left_lick_onset_idx)+1) : (max(left_lick_onset_idx)+resp_nBases);
reward_onset_idx = (max(right_lick_onset_idx)+1) : (max(right_lick_onset_idx)+reward_nBases);
leftforepaw_movement_idx = (max(reward_onset_idx)+1) : (max(reward_onset_idx)+movement_nBases);
rightforepaw_movement_idx = (max(leftforepaw_movement_idx)+1) : (max(leftforepaw_movement_idx)+movement_nBases);

predictor_non_movement_related_idx = 1:max(reward_onset_idx); % Non-movement related task predictors
predictor_both_movement_related_idx = [leftforepaw_movement_idx,rightforepaw_movement_idx]; % Non-movement related predictors
predictor_delay_right_lick_idx = [delay_left_middle_idx,delay_right_middle_idx,right_lick_onset_idx];

cell_coupling_idx = (max(rightforepaw_movement_idx)+1) : size(coeff{animal_num}{session_num},2);
assert(length(cell_coupling_idx) == size(coeff{animal_num}{session_num},2) - max(rightforepaw_movement_idx))

% Predictor groups for explained_deviance and p_value
task_tact_stim_onset_left_idx = 1;
task_tact_stim_onset_right_idx = 2;
task_delay_middle_left_idx = 3;
task_delay_middle_right_idx = 4;

% task_go_cue_idx = 5; % Removed gocue
task_left_lick_onset_idx = 5;
task_right_lick_onset_idx = 6;
task_reward_onset_idx = 7;
task_leftforepaw_movement_idx = 8;
task_rightforepaw_movement_idx = 9;

% Check if number of task grouped predictors is correct
assert(task_rightforepaw_movement_idx == (size(explained_deviance_all_reduction,2) - size(explained_deviance_all_reduction,1)),'Number of task grouped predictors is wrong')
last_idx = size(explained_deviance_all_reduction,2) - size(explained_deviance_all_reduction,1); % Demarcate the last task-variable predictor group index.

task_non_movement_related_idx = [1:task_reward_onset_idx]; % Non-movement related predictors
task_both_movement_related_idx = [task_leftforepaw_movement_idx,task_rightforepaw_movement_idx]; % Movement related predictors
task_delay_right_lick_idx = [ task_delay_middle_right_idx,task_right_lick_onset_idx];

task_cell_coupling_idx = (last_idx+1):size(explained_deviance_all_reduction,2);

% Test if number of task grouped predictors is correct
assert(size(task_cell_coupling_idx,2) == size(p_value,2))
assert(size(p_value_task1_var,2) == size(p_value_task1_var,1) +task_rightforepaw_movement_idx)

%% Setting up conditions to locate significant coupling cells
% Split up explained deviance matrix,p_value and coeff to differentiate between cell_coupling and non_cell coupling predictors
non_coupling_explained_deviance_all_reduction = explained_deviance_all_reduction(:,1:last_idx);
coupling_explained_deviance_all_reduction = explained_deviance_all_reduction(:,task_cell_coupling_idx);

non_coupling_p_value = p_value_task1_var(:,1:last_idx);
coupling_p_value = p_value_task1_var(:,task_cell_coupling_idx);

non_coupling_coeff = coeff{animal_num}{session_num}(:,1:max(rightforepaw_movement_idx));
coupling_coeff = coeff{animal_num}{session_num}(:,cell_coupling_idx);

assert( (size(non_coupling_coeff,2)+size(coupling_coeff,2)) == size(coeff{animal_num}{session_num},2),'Error in size of coupling coeff');
assert(size(coupling_coeff,1)*2 == size(coupling_coeff,2),'Error in size of coupling coeff');

%% Save useful variables into sessions for future analysis
all_resp_var_trial_concat{animal_num}{session_num} = resp_var_trial_concat; % Actual cell response
all_trial_type{animal_num}{session_num} = trial_type;
all_region_labels{animal_num}{session_num} = region_labels;
all_sessiondata{animal_num}{session_num} = SessionData;
all_GLM{animal_num}{session_num} = GLM_result;
all_p_value{animal_num}{session_num} = p_value;

all_coupling_explained_deviance_reduction{animal_num}{session_num} = coupling_explained_deviance_all_reduction;
all_coupling_p_value{animal_num}{session_num} = coupling_p_value;
all_coupling_coeff{animal_num}{session_num} = coupling_coeff;
all_coeff{animal_num}{session_num} = GLM_result.coeff;

assert(size(p_value,2) == size(all_coeff{animal_num}{session_num},1), 'Number of cells do not match - Possible incompletion of GLM or shuffling' )
all_ex_var{animal_num}{session_num} = GLM_result.explained_deviance_test;

all_predictor_labels{animal_num}{session_num} = predictor_labels;
%% Task predictors used to reconstruct activity
all_norm_tact_stim_onset_left_convolved_trial{animal_num}{session_num} =  norm_tact_stim_onset_left_convolved_trial;
all_norm_tact_stim_onset_right_convolved_trial{animal_num}{session_num} =  norm_tact_stim_onset_right_convolved_trial;
all_norm_delay_middle_left_convolved_trial{animal_num}{session_num} =  norm_delay_left_middle_convolved_trial;
all_norm_delay_middle_right_convolved_trial{animal_num}{session_num} =  norm_delay_right_middle_convolved_trial;
all_norm_lick_left_onset_convolved_trial{animal_num}{session_num} =  norm_lick_left_onset_convolved_trial;
all_norm_lick_right_onset_convolved_trial{animal_num}{session_num} =  norm_lick_right_onset_convolved_trial;
all_norm_reward_onset_convolved_trial{animal_num}{session_num} =  norm_reward_onset_convolved_trial;
all_norm_leftforepaw_movement_convolved_trial{animal_num}{session_num} =  norm_leftforepaw_movement_convolved_trial;
all_norm_rightforepaw_movement_convolved_trial{animal_num}{session_num} =  norm_rightforepaw_movement_convolved_trial;
all_norm_resp_task1_var_t1_trial{animal_num}{session_num} = norm_resp_task1_var_t1_trial;
all_norm_resp_task1_var_t2_trial{animal_num}{session_num} = norm_resp_task1_var_t2_trial;
all_task1_var_coupling_train{animal_num}{session_num} = task1_var_coupling_train;
all_task1_var_coupling_test{animal_num}{session_num} = task1_var_coupling_test;
all_resp_task1_var_test{animal_num}{session_num} = resp_task1_var_test;
all_resp_task1_var_train{animal_num}{session_num} = resp_task1_var_train;

% Valid cells are defined with explained deviance p value = 0 and coeff >0
all_valid_cells{animal_num}{session_num} = intersect(find(p_value == 0) , find(sum((GLM_result.coeff>0),2)>0)');

% Get trial CHOICE - labels for training and test data
% Note that we will use trial_types (1,2,3,4) but they will give different information. i.e. trial type 1 & 4 - Left, 2 % 3 - Right
all_trial_choice{animal_num}{session_num} = nan(length(all_sessiondata{animal_num}{session_num}.TrialTypes),1);
all_trial_choice{animal_num}{session_num}(find(trial_type == 1 | trial_type == 4)) = 1; % Left choice
all_trial_choice{animal_num}{session_num}(find(trial_type == 2 | trial_type == 3)) = 2; % Right choice

% To find significant predictor idx
for cell_num = 1:size(p_value,2)
    
    p_value_threshold = 0.05; % Or the critical value used to obtain the new threshold
    
    % Get sig cells for p value using FDR
    [sorted_coupling_p_value,b,sort_rank] = unique(coupling_p_value(cell_num,:));
    ranked_p_value_threshold = p_value_threshold * ([1:size(sorted_coupling_p_value,2)] / size(sorted_coupling_p_value,2));
    concat_ranked_p_value_threshold = [];
    for j = 1:length(sort_rank)
        concat_ranked_p_value_threshold = [concat_ranked_p_value_threshold ranked_p_value_threshold(sort_rank(j))];
    end
    coupling_p_value_sig_idx{cell_num} = find(coupling_p_value(cell_num,:) < concat_ranked_p_value_threshold);
    clear  ranked_p_value_threshold p_value_idx sorted_ranked_p_value_threshold concat_ranked_p_value_threshold ranked_p_value_threshold sorted_coupling_p_value
    
    [sorted_coupling_p_value,b,sort_rank] = unique(non_coupling_p_value(cell_num,:));
    ranked_p_value_threshold = p_value_threshold * ([1:size(sorted_coupling_p_value,2)] / size(sorted_coupling_p_value,2));
    concat_ranked_p_value_threshold = [];
    for j = 1:length(sort_rank)
        concat_ranked_p_value_threshold = [concat_ranked_p_value_threshold ranked_p_value_threshold(sort_rank(j))];
    end
    
    non_coupling_p_value_sig_idx{cell_num} = find(non_coupling_p_value(cell_num,:) < concat_ranked_p_value_threshold);
    clear  ranked_p_value_threshold p_value_idx sorted_ranked_p_value_threshold concat_ranked_p_value_threshold ranked_p_value_threshold sorted_coupling_p_value
    non_coupling_ev_reduction_sig_idx{cell_num} = find(non_coupling_explained_deviance_all_reduction(cell_num,:) < 0);
    non_coupling_coeff_sig_idx{cell_num} = find(non_coupling_coeff(cell_num,:) > 0); % For each predicted cell, which coeff is non-zero, i.e. cells that can be predicted. If coeff is = 0 or <0 means that behaviour unable to predict well.
    coupling_ev_reduction_sig_idx{cell_num} = find(coupling_explained_deviance_all_reduction(cell_num,:) < 0); % For each predicted cells, which coupling predictors EV reduction are <0
    
    coupling_coeff_sig_idx{cell_num} = []; % To find out which predictor cell is sig. to predict the predicted cell.
    
    % For +vely correlated cells - with positive coeff
    for i = 1:size(coupling_coeff(cell_num,:),2) % For each predictor, get the significant cell index
        predictor_idx = cell_coupling(i);
        % Append the cell_coupling_idx
        if coupling_coeff(cell_num,i) > 0 % If predictor cell is >0 , it means that it can be sig. used to predict predicted cell.
            coupling_coeff_sig_idx{cell_num} = [coupling_coeff_sig_idx{cell_num}, predictor_idx];
        end
    end
    
    coupling_coeff_sig_idx{cell_num} = unique(coupling_coeff_sig_idx{cell_num}); % Remove the duplicate sig cell idx
    neg_coupling_coeff_sig_idx{cell_num} = []; % To find out which predictor cell is sig. to predict the predicted cell.
end

% Get intersection of significant/valid: EV, P_value and Coeff
% Get significant cell indices
tact_left_sig_cells{animal_num}{session_num} = [];
tact_right_sig_cells{animal_num}{session_num} = [];
delay_middle_left_sig_cells{animal_num}{session_num}  = [];
delay_middle_right_sig_cells{animal_num}{session_num}  = [];
left_lick_sig_cells{animal_num}{session_num} = [];
right_lick_sig_cells{animal_num}{session_num} = [];

% Note all cells with pvalue > 0 will not be considered. i.e. they are taken as non-task relevant cell, i.e. the GLM is not applicable for the cell
for cell_num = 1:size(p_value,2)
    
    % Get significant and valid coupling predictors for each cell
    if p_value(cell_num) == 0 % To ensure that p_value for all task variables is significant
        % Get intersection of valid/sig p_value, coeff and explained_deviance
        coupling_sig_cells{animal_num}{session_num}{cell_num} = intersect(coupling_p_value_sig_idx{cell_num},coupling_coeff_sig_idx{cell_num});
    else coupling_sig_cells{animal_num}{session_num}{cell_num} = nan;
    end
    
    % Get indices of cells that are significant for specific behaviour/task variables
    % Tact stim left
    if p_value(cell_num) == 0 && (sum(non_coupling_ev_reduction_sig_idx{cell_num} == task_tact_stim_onset_left_idx) > 0) && (sum(non_coupling_p_value_sig_idx{cell_num} == task_tact_stim_onset_left_idx) > 0) && (length(intersect(non_coupling_coeff_sig_idx{cell_num},tact_stim_onset_left_idx))>0)
        tact_left_sig_cells{animal_num}{session_num} = [tact_left_sig_cells{animal_num}{session_num} , cell_num];
    end
    
    % Tact stim right
    if p_value(cell_num) == 0 && (sum(non_coupling_ev_reduction_sig_idx{cell_num} == task_tact_stim_onset_right_idx) > 0) && (sum(non_coupling_p_value_sig_idx{cell_num} == task_tact_stim_onset_right_idx) > 0) && (length(intersect(non_coupling_coeff_sig_idx{cell_num},tact_stim_onset_right_idx))>0)
        tact_right_sig_cells{animal_num}{session_num} = [tact_right_sig_cells{animal_num}{session_num} , cell_num];
    end
    
    % Delay left
    if p_value(cell_num) == 0 && (sum(non_coupling_ev_reduction_sig_idx{cell_num} == task_delay_middle_left_idx) > 0) && (sum(non_coupling_p_value_sig_idx{cell_num} == task_delay_middle_left_idx) > 0) && (length(intersect(non_coupling_coeff_sig_idx{cell_num},delay_left_middle_idx))>0)
        delay_middle_left_sig_cells{animal_num}{session_num} = [delay_middle_left_sig_cells{animal_num}{session_num} , cell_num];
    end
    
    % Delay right
    if p_value(cell_num) == 0 && (sum(non_coupling_ev_reduction_sig_idx{cell_num} == task_delay_middle_right_idx) > 0) && (sum(non_coupling_p_value_sig_idx{cell_num} == task_delay_middle_right_idx) > 0) && (length(intersect(non_coupling_coeff_sig_idx{cell_num},delay_right_middle_idx))>0)
        delay_middle_right_sig_cells{animal_num}{session_num} = [delay_middle_right_sig_cells{animal_num}{session_num} , cell_num];
    end
    
    % Left lick onset
    if p_value(cell_num) == 0 && (sum(non_coupling_ev_reduction_sig_idx{cell_num} == task_left_lick_onset_idx) > 0) && (sum(non_coupling_p_value_sig_idx{cell_num} == task_left_lick_onset_idx) > 0) && (length(intersect(non_coupling_coeff_sig_idx{cell_num},left_lick_onset_idx))>0)
        left_lick_sig_cells{animal_num}{session_num} = [left_lick_sig_cells{animal_num}{session_num} , cell_num];
    end
    % Right lick onset
    if p_value(cell_num) == 0 && (sum(non_coupling_ev_reduction_sig_idx{cell_num} == task_right_lick_onset_idx) > 0) && (sum(non_coupling_p_value_sig_idx{cell_num} == task_right_lick_onset_idx) > 0) && (length(intersect(non_coupling_coeff_sig_idx{cell_num},right_lick_onset_idx))>0)
        right_lick_sig_cells{animal_num}{session_num} = [right_lick_sig_cells{animal_num}{session_num} , cell_num];
    end
end

% Run some assertions to ensure that the GLM was done properly
assert(sum(GLM_result.explained_deviance_test ==0) <1,'Presence of Zeros in Full Ex Var. Model')

disp(['Total number of cells: ',num2str(size(p_value,2))])
disp(['Total number of valid cells: ',num2str(size(all_valid_cells{animal_num}{session_num},2))])
save(['GLM_concat_results.mat'],'-v7.3')