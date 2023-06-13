%% (Demo) Obtain full and partial models of the real and shuffled data.
% Load.
parameters.root_folder = uigetdir;
cd([parameters.root_folder,'/GLM_demo'])
load GLM_result.mat
load GLM_parameters.mat

fs_image = 9.352;
x_size = 500;
y_size = 500;
region = 8;
xy_pixel_ratio = 0.4;

% Full model.
B0 = GLM_result.B0;
coeff = GLM_result.coeff;
predict_train = GLM_result.predict_train;
predict_test = GLM_result.predict_test;
y_test = GLM_result.y_test;
y_hat_test = GLM_result.y_hat_test;
y_null_test = GLM_result.y_null_test;
L1_test = GLM_result.L1_test;
L0_test = GLM_result.L0_test;
LS_test = GLM_result.LS_test;
explained_deviance_test = GLM_result.explained_deviance_test;

% Make sure the GLM finished running
assert(length(B0) == size(resp_task1_var_test,2),'GLM not complete')

% Edited to remove go-cue
tact_stim_onset_left_idx = 1:stim_nBases;
tact_stim_onset_right_idx = (max(tact_stim_onset_left_idx)+1) : (max(tact_stim_onset_left_idx)+stim_nBases);
delay_left_middle_idx = (max(tact_stim_onset_right_idx)+1) : (max(tact_stim_onset_right_idx)+delay_nBases);
delay_right_middle_idx = (max(delay_left_middle_idx)+1) : (max(delay_left_middle_idx)+delay_nBases);
left_lick_onset_idx = (max(delay_right_middle_idx)+1) : (max(delay_right_middle_idx)+resp_nBases);
right_lick_onset_idx = (max(left_lick_onset_idx)+1) : (max(left_lick_onset_idx)+resp_nBases);
reward_onset_idx = (max(right_lick_onset_idx)+1) : (max(right_lick_onset_idx)+reward_nBases);
leftforepaw_movement_idx = (max(reward_onset_idx)+1) : (max(reward_onset_idx)+movement_nBases);
rightforepaw_movement_idx = (max(leftforepaw_movement_idx)+1) : (max(leftforepaw_movement_idx)+movement_nBases);

% To look at individual task predictors in the delay period (Note that this is the reverse)
delay_left_middle_1_idx = delay_left_middle_idx(1);
delay_left_middle_2_idx = delay_left_middle_idx(2);
delay_left_middle_3_idx = delay_left_middle_idx(3);
delay_left_middle_4_idx = delay_left_middle_idx(4);
delay_left_middle_5_idx = delay_left_middle_idx(5);

delay_right_middle_1_idx = delay_right_middle_idx(1);
delay_right_middle_2_idx = delay_right_middle_idx(2);
delay_right_middle_3_idx = delay_right_middle_idx(3);
delay_right_middle_4_idx = delay_right_middle_idx(4);
delay_right_middle_5_idx = delay_right_middle_idx(5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Concatenate all trials for task variables
task1_idx = union(task1_train_idx,task1_test_idx);
for trial_num = 1:length(task1_idx)
    task1_var_temp{trial_num} = [
        norm_tact_stim_onset_left_convolved_trial{task1_idx(trial_num)},...
        norm_tact_stim_onset_right_convolved_trial{task1_idx(trial_num)},...
        norm_delay_left_middle_convolved_trial{task1_idx(trial_num)},...
        norm_delay_right_middle_convolved_trial{task1_idx(trial_num)},...
        norm_lick_left_onset_convolved_trial{task1_idx(trial_num)},...
        norm_lick_right_onset_convolved_trial{task1_idx(trial_num)},...
        norm_reward_onset_convolved_trial{task1_idx(trial_num)},...
        norm_leftforepaw_movement_convolved_trial{task1_idx(trial_num)},...
        norm_rightforepaw_movement_convolved_trial{task1_idx(trial_num)}];
end
assert(size(task1_var_temp{1},2)==max(rightforepaw_movement_idx))

% Concatenate.
task1_var = [];
resp_task1_var = [];
for trial_num = 1:length(task1_idx)
    task1_var = [task1_var;task1_var_temp{trial_num}]; % Concatenate across the entire session
    resp_task1_var = [resp_task1_var;resp_var_trial_concat{task1_idx(trial_num)}];
end

% Add cell coupling predictors to task1_var
% Arrange t1 and t2 predictors for the same cell side by side
for cell_num = 1 : size(norm_resp_task1_var_t2,2)
    task1_var = [task1_var norm_resp_task1_var_t1(:,cell_num) norm_resp_task1_var_t2(:,cell_num)]; % Note that this is across the entire trial
end

% Get train and test frames.
task1_train_frame_all = [];
task1_test_frame_all = [];
for trial_num = 1:length(task1_train_idx)
    task1_train_frame_all = [task1_train_frame_all,(trial_begin_img(task1_train_idx(trial_num)) - 4.*round(fs_image)):trial_end_img(task1_train_idx(trial_num))];
end
for trial_num = 1:length(task1_test_idx)
    task1_test_frame_all = [task1_test_frame_all,(trial_begin_img(task1_test_idx(trial_num)) - 4.*round(fs_image)):trial_end_img(task1_test_idx(trial_num))];
end

assert(size(task1_train_frame_all,2)+size(task1_test_frame_all,2) == size(task1_var,1));

% To get the frame index from the actual time series
task1_train_test_frame_temp = zeros(1,trial_end_img(end));
task1_train_test_frame_temp(task1_train_frame_all) = 1;
task1_train_test_frame_temp(task1_test_frame_all) = 2;
task1_train_test_frame = task1_train_test_frame_temp(task1_train_test_frame_temp > 0);
train_frame = find(task1_train_test_frame == 1);
test_frame = find(task1_train_test_frame == 2);

%% Initialising variables to prelocate for speed
for hide = 1:1 % To hide chunk of code
    % Mean y_hat variables
    mean_y_hat_tact_stim_onset_left = nan(length(GLM_result.B0),1);
    mean_y_hat_tact_stim_onset_right = nan(length(GLM_result.B0),1);
    mean_y_hat_delay_middle_left = nan(length(GLM_result.B0),1);
    mean_y_hat_delay_middle_right = nan(length(GLM_result.B0),1);
    % mean_y_hat_go_cue = nan(length(GLM_result.B0),1);
    mean_y_hat_left_lick_onset = nan(length(GLM_result.B0),1);
    mean_y_hat_right_lick_onset = nan(length(GLM_result.B0),1);
    mean_y_hat_reward_onset = nan(length(GLM_result.B0),1);
    mean_y_hat_leftforepaw_movement = nan(length(GLM_result.B0),1);
    mean_y_hat_rightforepaw_movement = nan(length(GLM_result.B0),1);
    
    mean_y_hat_delay_middle_left_1 = nan(length(GLM_result.B0),1);
    mean_y_hat_delay_middle_left_2 = nan(length(GLM_result.B0),1);
    mean_y_hat_delay_middle_left_3 = nan(length(GLM_result.B0),1);
    mean_y_hat_delay_middle_left_4 = nan(length(GLM_result.B0),1);
    mean_y_hat_delay_middle_right_1 = nan(length(GLM_result.B0),1);
    mean_y_hat_delay_middle_right_2 = nan(length(GLM_result.B0),1);
    mean_y_hat_delay_middle_right_3 = nan(length(GLM_result.B0),1);
    mean_y_hat_delay_middle_right_4 = nan(length(GLM_result.B0),1);
    
    % y_hat_variables
    y_hat_tact_stim_onset_left_removed = nan(length(GLM_result.B0),size(task1_var,1));
    y_hat_tact_stim_onset_left_removed_test= nan(length(GLM_result.B0),size(test_frame,2));
    L1_tact_stim_onset_left_removed_test = nan(length(GLM_result.B0),1);
    explained_deviance_tact_stim_onset_left_removed_test = nan(length(GLM_result.B0),1);
    
    y_hat_tact_stim_onset_right_removed = nan(length(GLM_result.B0),size(task1_var,1));
    y_hat_tact_stim_onset_right_removed_test= nan(length(GLM_result.B0),size(test_frame,2));
    L1_tact_stim_onset_right_removed_test = nan(length(GLM_result.B0),1);
    explained_deviance_tact_stim_onset_right_removed_test = nan(length(GLM_result.B0),1);
    
    y_hat_delay_middle_left_removed = nan(length(GLM_result.B0),size(task1_var,1));
    y_hat_delay_middle_left_removed_test= nan(length(GLM_result.B0),size(test_frame,2));
    L1_delay_middle_left_removed_test = nan(length(GLM_result.B0),1);
    explained_deviance_delay_middle_left_removed_test = nan(length(GLM_result.B0),1);
    
    y_hat_delay_middle_right_removed = nan(length(GLM_result.B0),size(task1_var,1));
    y_hat_delay_middle_right_removed_test= nan(length(GLM_result.B0),size(test_frame,2));
    L1_delay_middle_right_removed_test = nan(length(GLM_result.B0),1);
    explained_deviance_delay_middle_right_removed_test = nan(length(GLM_result.B0),1);
    %
    % y_hat_go_cue_removed = nan(length(GLM_result.B0),size(task1_var,1));
    % y_hat_go_cue_removed_test= nan(length(GLM_result.B0),size(test_frame,2));
    % L1_go_cue_removed_test = nan(length(GLM_result.B0),1);
    % explained_deviance_go_cue_removed_test = nan(length(GLM_result.B0),1);
    
    y_hat_left_lick_onset_removed = nan(length(GLM_result.B0),size(task1_var,1));
    y_hat_left_lick_onset_removed_test= nan(length(GLM_result.B0),size(test_frame,2));
    L1_left_lick_onset_removed_test = nan(length(GLM_result.B0),1);
    explained_deviance_left_lick_onset_removed_test = nan(length(GLM_result.B0),1);
    
    y_hat_right_lick_onset_removed = nan(length(GLM_result.B0),size(task1_var,1));
    y_hat_right_lick_onset_removed_test= nan(length(GLM_result.B0),size(test_frame,2));
    L1_right_lick_onset_removed_test = nan(length(GLM_result.B0),1);
    explained_deviance_right_lick_onset_removed_test = nan(length(GLM_result.B0),1);
    
    y_hat_reward_onset_removed = nan(length(GLM_result.B0),size(task1_var,1));
    y_hat_reward_onset_removed_test= nan(length(GLM_result.B0),size(test_frame,2));
    L1_reward_onset_removed_test = nan(length(GLM_result.B0),1);
    explained_deviance_reward_onset_removed_test = nan(length(GLM_result.B0),1);
    
    y_hat_leftforepaw_movement_removed = nan(length(GLM_result.B0),size(task1_var,1));
    y_hat_leftforepaw_movement_removed_test= nan(length(GLM_result.B0),size(test_frame,2));
    L1_leftforepaw_movement_removed_test = nan(length(GLM_result.B0),1);
    explained_deviance_leftforepaw_movement_removed_test = nan(length(GLM_result.B0),1);
    
    y_hat_rightforepaw_movement_removed = nan(length(GLM_result.B0),size(task1_var,1));
    y_hat_rightforepaw_movement_removed_test = nan(length(GLM_result.B0),size(test_frame,2));
    L1_rightforepaw_movement_removed_test = nan(length(GLM_result.B0),1);
    explained_deviance_rightforepaw_movement_removed_test = nan(length(GLM_result.B0),1);
    
    % y_hat_delay_middle_left_1_removed = nan(length(GLM_result.B0),size(task1_var,1));
    % y_hat_delay_middle_left_1_removed_test = nan(length(GLM_result.B0),size(test_frame,2));
    % L1_delay_middle_left_1_removed_test = nan(length(GLM_result.B0),1);
    % explained_deviance_delay_middle_left_1_removed_test = nan(length(GLM_result.B0),1);
    %
    % y_hat_delay_middle_left_2_removed = nan(length(GLM_result.B0),size(task1_var,1));
    % y_hat_delay_middle_left_2_removed_test = nan(length(GLM_result.B0),size(test_frame,2));
    % L1_delay_middle_left_2_removed_test = nan(length(GLM_result.B0),1);
    % explained_deviance_delay_middle_left_2_removed_test = nan(length(GLM_result.B0),1);
    %
    % y_hat_delay_middle_left_3_removed = nan(length(GLM_result.B0),size(task1_var,1));
    % y_hat_delay_middle_left_3_removed_test = nan(length(GLM_result.B0),size(test_frame,2));
    % L1_delay_middle_left_3_removed_test = nan(length(GLM_result.B0),1);
    % explained_deviance_delay_middle_left_3_removed_test = nan(length(GLM_result.B0),1);
    %
    % y_hat_delay_middle_left_4_removed = nan(length(GLM_result.B0),size(task1_var,1));
    % y_hat_delay_middle_left_4_removed_test = nan(length(GLM_result.B0),size(test_frame,2));
    % L1_delay_middle_left_4_removed_test = nan(length(GLM_result.B0),1);
    % explained_deviance_delay_middle_left_4_removed_test = nan(length(GLM_result.B0),1);
    %
    % y_hat_delay_middle_right_1_removed = nan(length(GLM_result.B0),size(task1_var,1));
    % y_hat_delay_middle_right_1_removed_test = nan(length(GLM_result.B0),size(test_frame,2));
    % L1_delay_middle_right_1_removed_test = nan(length(GLM_result.B0),1);
    % explained_deviance_delay_middle_right_1_removed_test = nan(length(GLM_result.B0),1);
    %
    % y_hat_delay_middle_right_2_removed = nan(length(GLM_result.B0),size(task1_var,1));
    % y_hat_delay_middle_right_2_removed_test = nan(length(GLM_result.B0),size(test_frame,2));
    % L1_delay_middle_right_2_removed_test = nan(length(GLM_result.B0),1);
    % explained_deviance_delay_middle_right_2_removed_test = nan(length(GLM_result.B0),1);
    %
    % y_hat_delay_middle_right_3_removed = nan(length(GLM_result.B0),size(task1_var,1));
    % y_hat_delay_middle_right_3_removed_test = nan(length(GLM_result.B0),size(test_frame,2));
    % L1_delay_middle_right_3_removed_test = nan(length(GLM_result.B0),1);
    % explained_deviance_delay_middle_right_3_removed_test = nan(length(GLM_result.B0),1);
    %
    % y_hat_delay_middle_right_right_4_removed = nan(length(GLM_result.B0),size(task1_var,1));
    % y_hat_delay_middle_right_4_removed_test = nan(length(GLM_result.B0),size(test_frame,2));
    % L1_delay_middle_right_4_removed_test = nan(length(GLM_result.B0),1);
    % explained_deviance_delay_middle_right_4_removed_test = nan(length(GLM_result.B0),1);
    
    % y_hat_delay_middle_5_removed = nan(length(GLM_result.B0),size(task1_var,1));
    % y_hat_delay_middle_5_removed_test = nan(length(GLM_result.B0),size(test_frame,2));
    % L1_delay_middle_5_removed_test = nan(length(GLM_result.B0),1);
    % explained_deviance_delay_middle_5_removed_test = nan(length(GLM_result.B0),1);
    %
    % y_hat_delay_middle_6_removed = nan(length(GLM_result.B0),size(task1_var,1));
    % y_hat_delay_middle_6_removed_test = nan(length(GLM_result.B0),size(test_frame,2));
    % L1_delay_middle_6_removed_test = nan(length(GLM_result.B0),1);
    % explained_deviance_delay_middle_6_removed_test = nan(length(GLM_result.B0),1);
    %
    % y_hat_delay_middle_7_removed = nan(length(GLM_result.B0),size(task1_var,1));
    % y_hat_delay_middle_7_removed_test = nan(length(GLM_result.B0),size(test_frame,2));
    % L1_delay_middle_7_removed_test = nan(length(GLM_result.B0),1);
    % explained_deviance_delay_middle_7_removed_test = nan(length(GLM_result.B0),1);
    %
    % y_hat_delay_middle_8_removed = nan(length(GLM_result.B0),size(task1_var,1));
    % y_hat_delay_middle_8_removed_test = nan(length(GLM_result.B0),size(test_frame,2));
    % L1_delay_middle_8_removed_test = nan(length(GLM_result.B0),1);
    % explained_deviance_delay_middle_8_removed_test = nan(length(GLM_result.B0),1);
    
    
    % Initialise cell coupling variables using number of cells. Note that T1 and T2 coupling will be removed as a pair
    mean_y_hat_cell_coupling = nan(length(GLM_result.B0),length(GLM_result.B0)); % (Cell coupling predictors x predicted cells x 1)
    y_hat_cell_coupling_removed = nan(length(GLM_result.B0),length(GLM_result.B0),size(task1_var,1)); % (Cell coupling predictors x predicted cells x bins)
    y_hat_cell_coupling_removed_test = nan(length(GLM_result.B0),length(GLM_result.B0),size(test_frame,2));
    L1_cell_coupling_removed_test = nan(length(GLM_result.B0),length(GLM_result.B0),1);
    explained_deviance_cell_coupling_removed_test = nan(length(GLM_result.B0),length(GLM_result.B0),1);
    
    % Get mean of cell coupling predictors using for loop
    for predictor_cell_num = 1 : max(unique(predictor_labels(~isnan(predictor_labels))))
        % Get predictor labels for t1 and t2
        predictor_idx = find(predictor_labels(~isnan(predictor_labels)) == predictor_cell_num) + max(find(isnan(predictor_labels)));
        % I.e. Taking the mean of the predicted activity at each predictor cell index
        mean_y_hat_cell_coupling(predictor_cell_num,:) = mean(exp(task1_var(:,predictor_idx)*coeff(:,predictor_idx)'));
    end
    
    % Scaling factor marginalizing out the effect of the other variables.
    % Getting the mean of the predicted activity
    mean_y_hat_tact_stim_onset_left = mean(exp(task1_var(:,tact_stim_onset_left_idx)*coeff(:,tact_stim_onset_left_idx)'));
    mean_y_hat_tact_stim_onset_right = mean(exp(task1_var(:,tact_stim_onset_right_idx)*coeff(:,tact_stim_onset_right_idx)'));
    mean_y_hat_delay_middle_left = mean(exp(task1_var(:,delay_left_middle_idx)*coeff(:,delay_left_middle_idx)'));
    mean_y_hat_delay_middle_right = mean(exp(task1_var(:,delay_right_middle_idx)*coeff(:,delay_right_middle_idx)'));
    
    % mean_y_hat_go_cue = mean(exp(task1_var(:,go_cue_idx)*coeff(:,go_cue_idx)'));
    mean_y_hat_left_lick_onset = mean(exp(task1_var(:,left_lick_onset_idx)*coeff(:,left_lick_onset_idx)'));
    mean_y_hat_right_lick_onset = mean(exp(task1_var(:,right_lick_onset_idx)*coeff(:,right_lick_onset_idx)'));
    mean_y_hat_reward_onset = mean(exp(task1_var(:,reward_onset_idx)*coeff(:,reward_onset_idx)'));
    mean_y_hat_leftforepaw_movement = mean(exp(task1_var(:,leftforepaw_movement_idx)*coeff(:,leftforepaw_movement_idx)'));
    mean_y_hat_rightforepaw_movement = mean(exp(task1_var(:,rightforepaw_movement_idx)*coeff(:,rightforepaw_movement_idx)'));
    %
    % mean_y_hat_delay_middle_left_1 = mean(exp(task1_var(:,delay_left_middle_1_idx)*coeff(:,delay_left_middle_1_idx)'));
    % mean_y_hat_delay_middle_left_2 = mean(exp(task1_var(:,delay_left_middle_2_idx)*coeff(:,delay_left_middle_2_idx)'));
    % mean_y_hat_delay_middle_left_3 = mean(exp(task1_var(:,delay_left_middle_3_idx)*coeff(:,delay_left_middle_3_idx)'));
    % mean_y_hat_delay_middle_left_4 = mean(exp(task1_var(:,delay_left_middle_4_idx)*coeff(:,delay_left_middle_4_idx)'));
    % mean_y_hat_delay_middle_right_1 = mean(exp(task1_var(:,delay_right_middle_1_idx)*coeff(:,delay_right_middle_1_idx)'));
    % mean_y_hat_delay_middle_right_2 = mean(exp(task1_var(:,delay_right_middle_2_idx)*coeff(:,delay_right_middle_2_idx)'));
    % mean_y_hat_delay_middle_right_3 = mean(exp(task1_var(:,delay_right_middle_3_idx)*coeff(:,delay_right_middle_3_idx)'));
    % mean_y_hat_delay_middle_right_4 = mean(exp(task1_var(:,delay_right_middle_4_idx)*coeff(:,delay_right_middle_4_idx)'));
    % mean_y_hat_delay_middle_5 = mean(exp(task1_var(:,delay_middle_5_idx)*coeff(:,delay_middle_5_idx)'));
    % mean_y_hat_delay_middle_6 = mean(exp(task1_var(:,delay_middle_6_idx)*coeff(:,delay_middle_6_idx)'));
    % mean_y_hat_delay_middle_7 = mean(exp(task1_var(:,delay_middle_7_idx)*coeff(:,delay_middle_7_idx)'));
    % mean_y_hat_delay_middle_8 = mean(exp(task1_var(:,delay_middle_8_idx)*coeff(:,delay_middle_8_idx)'));
    clear y_test y_null_test L0_test LS_test
    for cell_num = 1: length(GLM_result.B0)
        % Model assessment by measuring explained deviance based on Benjamin et al 2018.
        y_test(cell_num,:) = resp_task1_var_test(:,cell_num);
        y_null_test(cell_num) = mean(resp_task1_var_test(:,cell_num),1);
        L0_test(cell_num) = sum(y_test(cell_num,:).*log(eps + y_null_test(cell_num)) - y_null_test(cell_num)); % Null model.
        LS_test(cell_num) = sum(y_test(cell_num,:).*log(eps + y_test(cell_num,:)) - y_test(cell_num,:)); % Saturated model.
        
    end
end

%% Get predicted activity after removing the specific predictor
assert( max(find(isnan(predictor_labels))) == max(rightforepaw_movement_idx));
start_coupling_idx = max(find(isnan(predictor_labels))) + 1;
end_coupling_idx = size(task1_var,2);

% Remove cell_coupling using (Predictor cells (Activity removed) x Predicted cells x Bins) size for y_hat_cell_coupling_removed
parfor predictor_cell_num = 1 : max(unique(predictor_labels(~isnan(predictor_labels))))
    predictor_idx = find(predictor_labels(~isnan(predictor_labels)) == predictor_cell_num) + max(find(isnan(predictor_labels)));
    % Get the index of all cell coupling predictors except the one being removed
    all_predictors_except = setdiff(start_coupling_idx:end_coupling_idx, predictor_idx)
    
    % Get 3D tensor of y_hat, with the first dim. being the predictor cell that is removed
    y_hat_cell_coupling_removed(predictor_cell_num,:,:) = ...
        (exp(task1_var(:,all_predictors_except)*coeff(:,all_predictors_except)').*...
        mean_y_hat_cell_coupling(predictor_cell_num,:).*...
        exp(task1_var(:,tact_stim_onset_left_idx)*coeff(:,tact_stim_onset_left_idx)').*...
        exp(task1_var(:,tact_stim_onset_right_idx)*coeff(:,tact_stim_onset_right_idx)').*...
        exp(task1_var(:,delay_left_middle_idx)*coeff(:,delay_left_middle_idx)').*...
        exp(task1_var(:,delay_right_middle_idx)*coeff(:,delay_right_middle_idx)').*...
        exp(task1_var(:,left_lick_onset_idx)*coeff(:,left_lick_onset_idx)').*...
        exp(task1_var(:,right_lick_onset_idx)*coeff(:,right_lick_onset_idx)').*...
        exp(task1_var(:,reward_onset_idx)*coeff(:,reward_onset_idx)').*...
        exp(task1_var(:,leftforepaw_movement_idx)*coeff(:,leftforepaw_movement_idx)').*...
        exp(task1_var(:,rightforepaw_movement_idx)*coeff(:,rightforepaw_movement_idx)').*...
        exp(B0))';
end

% Remove tactile stimulus onset left.
matrix_y_hat_tact_stim_onset_left_removed = ...
    mean_y_hat_tact_stim_onset_left.*...
    exp(task1_var(:,start_coupling_idx:end_coupling_idx)*coeff(:,start_coupling_idx:end_coupling_idx)').*...
    exp(task1_var(:,tact_stim_onset_right_idx)*coeff(:,tact_stim_onset_right_idx)').*...
    exp(task1_var(:,delay_left_middle_idx)*coeff(:,delay_left_middle_idx)').*...
    exp(task1_var(:,delay_right_middle_idx)*coeff(:,delay_right_middle_idx)').*...
    exp(task1_var(:,left_lick_onset_idx)*coeff(:,left_lick_onset_idx)').*...
    exp(task1_var(:,right_lick_onset_idx)*coeff(:,right_lick_onset_idx)').*...
    exp(task1_var(:,reward_onset_idx)*coeff(:,reward_onset_idx)').*...
    exp(task1_var(:,leftforepaw_movement_idx)*coeff(:,leftforepaw_movement_idx)').*...
    exp(task1_var(:,rightforepaw_movement_idx)*coeff(:,rightforepaw_movement_idx)').*...
    exp(B0);
y_hat_tact_stim_onset_left_removed = matrix_y_hat_tact_stim_onset_left_removed';

% Remove tactile stimulus onset right.
matrix_y_hat_tact_stim_onset_right_removed = ...
    exp(task1_var(:,start_coupling_idx:end_coupling_idx)*coeff(:,start_coupling_idx:end_coupling_idx)').*...
    exp(task1_var(:,tact_stim_onset_left_idx)*coeff(:,tact_stim_onset_left_idx)').*...
    mean_y_hat_tact_stim_onset_right.*...
    exp(task1_var(:,delay_left_middle_idx)*coeff(:,delay_left_middle_idx)').*...
    exp(task1_var(:,delay_right_middle_idx)*coeff(:,delay_right_middle_idx)').*...
    exp(task1_var(:,left_lick_onset_idx)*coeff(:,left_lick_onset_idx)').*...
    exp(task1_var(:,right_lick_onset_idx)*coeff(:,right_lick_onset_idx)').*...
    exp(task1_var(:,reward_onset_idx)*coeff(:,reward_onset_idx)').*...
    exp(task1_var(:,leftforepaw_movement_idx)*coeff(:,leftforepaw_movement_idx)').*...
    exp(task1_var(:,rightforepaw_movement_idx)*coeff(:,rightforepaw_movement_idx)').*...
    exp(B0);
y_hat_tact_stim_onset_right_removed = matrix_y_hat_tact_stim_onset_right_removed';

% Remove delay middle left.
matrix_y_hat_delay_middle_left_removed = ...
    exp(task1_var(:,start_coupling_idx:end_coupling_idx)*coeff(:,start_coupling_idx:end_coupling_idx)').*...
    exp(task1_var(:,tact_stim_onset_left_idx)*coeff(:,tact_stim_onset_left_idx)').*...
    exp(task1_var(:,tact_stim_onset_right_idx)*coeff(:,tact_stim_onset_right_idx)').*...
    mean_y_hat_delay_middle_left.*...
    exp(task1_var(:,delay_right_middle_idx)*coeff(:,delay_right_middle_idx)').*...
    exp(task1_var(:,left_lick_onset_idx)*coeff(:,left_lick_onset_idx)').*...
    exp(task1_var(:,right_lick_onset_idx)*coeff(:,right_lick_onset_idx)').*...
    exp(task1_var(:,reward_onset_idx)*coeff(:,reward_onset_idx)').*...
    exp(task1_var(:,leftforepaw_movement_idx)*coeff(:,leftforepaw_movement_idx)').*...
    exp(task1_var(:,rightforepaw_movement_idx)*coeff(:,rightforepaw_movement_idx)').*...
    exp(B0);
y_hat_delay_middle_left_removed = matrix_y_hat_delay_middle_left_removed';

% Remove delay middle right
matrix_y_hat_delay_middle_right_removed = ...
    exp(task1_var(:,start_coupling_idx:end_coupling_idx)*coeff(:,start_coupling_idx:end_coupling_idx)').*...
    exp(task1_var(:,tact_stim_onset_left_idx)*coeff(:,tact_stim_onset_left_idx)').*...
    exp(task1_var(:,tact_stim_onset_right_idx)*coeff(:,tact_stim_onset_right_idx)').*...
    exp(task1_var(:,delay_left_middle_idx)*coeff(:,delay_left_middle_idx)').*...
    mean_y_hat_delay_middle_right.*...
    exp(task1_var(:,left_lick_onset_idx)*coeff(:,left_lick_onset_idx)').*...
    exp(task1_var(:,right_lick_onset_idx)*coeff(:,right_lick_onset_idx)').*...
    exp(task1_var(:,reward_onset_idx)*coeff(:,reward_onset_idx)').*...
    exp(task1_var(:,leftforepaw_movement_idx)*coeff(:,leftforepaw_movement_idx)').*...
    exp(task1_var(:,rightforepaw_movement_idx)*coeff(:,rightforepaw_movement_idx)').*...
    exp(B0);
y_hat_delay_middle_right_removed = matrix_y_hat_delay_middle_right_removed';

% Remove left lick onset.
matrix_y_hat_left_lick_onset_removed = ...
    exp(task1_var(:,start_coupling_idx:end_coupling_idx)*coeff(:,start_coupling_idx:end_coupling_idx)').*...
    exp(task1_var(:,tact_stim_onset_left_idx)*coeff(:,tact_stim_onset_left_idx)').*...
    exp(task1_var(:,tact_stim_onset_right_idx)*coeff(:,tact_stim_onset_right_idx)').*...
    exp(task1_var(:,delay_left_middle_idx)*coeff(:,delay_left_middle_idx)').*...
    exp(task1_var(:,delay_right_middle_idx)*coeff(:,delay_right_middle_idx)').*...
    mean_y_hat_left_lick_onset.*...
    exp(task1_var(:,right_lick_onset_idx)*coeff(:,right_lick_onset_idx)').*...
    exp(task1_var(:,reward_onset_idx)*coeff(:,reward_onset_idx)').*...
    exp(task1_var(:,leftforepaw_movement_idx)*coeff(:,leftforepaw_movement_idx)').*...
    exp(task1_var(:,rightforepaw_movement_idx)*coeff(:,rightforepaw_movement_idx)').*...
    exp(B0);
y_hat_left_lick_onset_removed = matrix_y_hat_left_lick_onset_removed';

% Remove right lick onset.
matrix_y_hat_right_lick_onset_removed = ...
    exp(task1_var(:,start_coupling_idx:end_coupling_idx)*coeff(:,start_coupling_idx:end_coupling_idx)').*...
    exp(task1_var(:,tact_stim_onset_left_idx)*coeff(:,tact_stim_onset_left_idx)').*...
    exp(task1_var(:,tact_stim_onset_right_idx)*coeff(:,tact_stim_onset_right_idx)').*...
    exp(task1_var(:,delay_left_middle_idx)*coeff(:,delay_left_middle_idx)').*...
    exp(task1_var(:,delay_right_middle_idx)*coeff(:,delay_right_middle_idx)').*...
    exp(task1_var(:,left_lick_onset_idx)*coeff(:,left_lick_onset_idx)').*...
    mean_y_hat_right_lick_onset.*...
    exp(task1_var(:,reward_onset_idx)*coeff(:,reward_onset_idx)').*...
    exp(task1_var(:,leftforepaw_movement_idx)*coeff(:,leftforepaw_movement_idx)').*...
    exp(task1_var(:,rightforepaw_movement_idx)*coeff(:,rightforepaw_movement_idx)').*...
    exp(B0);
y_hat_right_lick_onset_removed = matrix_y_hat_right_lick_onset_removed';

% Remove reward onset.
matrix_y_hat_reward_onset_removed = ...
    exp(task1_var(:,start_coupling_idx:end_coupling_idx)*coeff(:,start_coupling_idx:end_coupling_idx)').*...
    exp(task1_var(:,tact_stim_onset_left_idx)*coeff(:,tact_stim_onset_left_idx)').*...
    exp(task1_var(:,tact_stim_onset_right_idx)*coeff(:,tact_stim_onset_right_idx)').*...
    exp(task1_var(:,delay_left_middle_idx)*coeff(:,delay_left_middle_idx)').*...
    exp(task1_var(:,delay_right_middle_idx)*coeff(:,delay_right_middle_idx)').*...
    exp(task1_var(:,left_lick_onset_idx)*coeff(:,left_lick_onset_idx)').*...
    exp(task1_var(:,right_lick_onset_idx)*coeff(:,right_lick_onset_idx)').*...
    mean_y_hat_reward_onset.*...
    exp(task1_var(:,leftforepaw_movement_idx)*coeff(:,leftforepaw_movement_idx)').*...
    exp(task1_var(:,rightforepaw_movement_idx)*coeff(:,rightforepaw_movement_idx)').*...
    exp(B0);
y_hat_reward_onset_removed = matrix_y_hat_reward_onset_removed';

% Remove left forepaw movement.
matrix_y_hat_leftforepaw_movement_removed = ...
    exp(task1_var(:,start_coupling_idx:end_coupling_idx)*coeff(:,start_coupling_idx:end_coupling_idx)').*...
    exp(task1_var(:,tact_stim_onset_left_idx)*coeff(:,tact_stim_onset_left_idx)').*...
    exp(task1_var(:,tact_stim_onset_right_idx)*coeff(:,tact_stim_onset_right_idx)').*...
    exp(task1_var(:,delay_left_middle_idx)*coeff(:,delay_left_middle_idx)').*...
    exp(task1_var(:,delay_right_middle_idx)*coeff(:,delay_right_middle_idx)').*...
    exp(task1_var(:,left_lick_onset_idx)*coeff(:,left_lick_onset_idx)').*...
    exp(task1_var(:,right_lick_onset_idx)*coeff(:,right_lick_onset_idx)').*...
    exp(task1_var(:,reward_onset_idx)*coeff(:,reward_onset_idx)').*...
    mean_y_hat_leftforepaw_movement.*...
    exp(task1_var(:,rightforepaw_movement_idx)*coeff(:,rightforepaw_movement_idx)').*...
    exp(B0);
y_hat_leftforepaw_movement_removed  = matrix_y_hat_leftforepaw_movement_removed';

% Remove right forepaw movement.
matrix_y_hat_rightforepaw_movement_removed = ...
    exp(task1_var(:,start_coupling_idx:end_coupling_idx)*coeff(:,start_coupling_idx:end_coupling_idx)').*...
    exp(task1_var(:,tact_stim_onset_left_idx)*coeff(:,tact_stim_onset_left_idx)').*...
    exp(task1_var(:,tact_stim_onset_right_idx)*coeff(:,tact_stim_onset_right_idx)').*...
    exp(task1_var(:,delay_left_middle_idx)*coeff(:,delay_left_middle_idx)').*...
    exp(task1_var(:,delay_right_middle_idx)*coeff(:,delay_right_middle_idx)').*...
    exp(task1_var(:,left_lick_onset_idx)*coeff(:,left_lick_onset_idx)').*...
    exp(task1_var(:,right_lick_onset_idx)*coeff(:,right_lick_onset_idx)').*...
    exp(task1_var(:,reward_onset_idx)*coeff(:,reward_onset_idx)').*...
    exp(task1_var(:,leftforepaw_movement_idx)*coeff(:,leftforepaw_movement_idx)').*...
    mean_y_hat_rightforepaw_movement.*...
    exp(B0);
y_hat_rightforepaw_movement_removed  = matrix_y_hat_rightforepaw_movement_removed';

%% Get explained deviance of removed model
clear explained_deviance_cell_coupling_removed_test
parfor predictor_cell_num = 1 : max(unique(predictor_labels(~isnan(predictor_labels))))
    y_hat_cell_coupling_removed_test = squeeze(y_hat_cell_coupling_removed(predictor_cell_num, :, test_frame)); % Only take into account the tests trials
    L1_cell_coupling_removed_test = sum(y_test.*log(eps + y_hat_cell_coupling_removed_test) - y_hat_cell_coupling_removed_test,2); % Partial model.
    
    % Arrange predictor cell_num removed idx in columns
    explained_deviance_cell_coupling_removed_test(predictor_cell_num,:) = 1 - (LS_test - L1_cell_coupling_removed_test')./(LS_test - L0_test);
end

y_hat_tact_stim_onset_left_removed_test = y_hat_tact_stim_onset_left_removed(:,test_frame); % Only take into account the tests trials
L1_tact_stim_onset_left_removed_test = sum(y_test.*log(eps + y_hat_tact_stim_onset_left_removed_test) - y_hat_tact_stim_onset_left_removed_test,2); % Partial model.
explained_deviance_tact_stim_onset_left_removed_test = 1 - (LS_test - L1_tact_stim_onset_left_removed_test')./(LS_test - L0_test);

y_hat_tact_stim_onset_right_removed_test = y_hat_tact_stim_onset_right_removed(:,test_frame);
L1_tact_stim_onset_right_removed_test = sum(y_test.*log(eps + y_hat_tact_stim_onset_right_removed_test) - y_hat_tact_stim_onset_right_removed_test,2); % Partial model.
explained_deviance_tact_stim_onset_right_removed_test = 1 - (LS_test - L1_tact_stim_onset_right_removed_test')./(LS_test - L0_test);

y_hat_delay_middle_left_removed_test = y_hat_delay_middle_left_removed(:,test_frame);
L1_delay_middle_left_removed_test = sum(y_test.*log(eps + y_hat_delay_middle_left_removed_test) - y_hat_delay_middle_left_removed_test,2); % Partial model.
explained_deviance_delay_middle_left_removed_test = 1 - (LS_test - L1_delay_middle_left_removed_test')./(LS_test - L0_test);

y_hat_delay_middle_right_removed_test = y_hat_delay_middle_right_removed(:,test_frame);
L1_delay_middle_right_removed_test = sum(y_test.*log(eps + y_hat_delay_middle_right_removed_test) - y_hat_delay_middle_right_removed_test,2); % Partial model.
explained_deviance_delay_middle_right_removed_test = 1 - (LS_test - L1_delay_middle_right_removed_test')./(LS_test - L0_test);

y_hat_left_lick_onset_removed_test = y_hat_left_lick_onset_removed(:,test_frame);
L1_left_lick_onset_removed_test = sum(y_test.*log(eps + y_hat_left_lick_onset_removed_test) - y_hat_left_lick_onset_removed_test,2); % Partial model.
explained_deviance_left_lick_onset_removed_test = 1 - (LS_test - L1_left_lick_onset_removed_test')./(LS_test - L0_test);

y_hat_right_lick_onset_removed_test = y_hat_right_lick_onset_removed(:,test_frame);
L1_right_lick_onset_removed_test = sum(y_test.*log(eps + y_hat_right_lick_onset_removed_test) - y_hat_right_lick_onset_removed_test,2); % Partial model.
explained_deviance_right_lick_onset_removed_test = 1 - (LS_test - L1_right_lick_onset_removed_test')./(LS_test - L0_test);

y_hat_reward_onset_removed_test = y_hat_reward_onset_removed(:,test_frame);
L1_reward_onset_removed_test = sum(y_test.*log(eps + y_hat_reward_onset_removed_test) - y_hat_reward_onset_removed_test,2); % Partial model.
explained_deviance_reward_onset_removed_test = 1 - (LS_test - L1_reward_onset_removed_test')./(LS_test - L0_test);

y_hat_leftforepaw_movement_removed_test = y_hat_leftforepaw_movement_removed(:,test_frame);
L1_leftforepaw_movement_removed_test = sum(y_test.*log(eps + y_hat_leftforepaw_movement_removed_test) - y_hat_leftforepaw_movement_removed_test,2); % Partial model.
explained_deviance_leftforepaw_movement_removed_test = 1 - (LS_test - L1_leftforepaw_movement_removed_test')./(LS_test - L0_test);

y_hat_rightforepaw_movement_removed_test = y_hat_rightforepaw_movement_removed(:,test_frame);
L1_rightforepaw_movement_removed_test = sum(y_test.*log(eps + y_hat_rightforepaw_movement_removed_test) - y_hat_rightforepaw_movement_removed_test,2); % Partial model.
explained_deviance_rightforepaw_movement_removed_test = 1 - (LS_test - L1_rightforepaw_movement_removed_test')./(LS_test - L0_test);

% Assess the partial models. (Cells x Number of removed predictor groups)
explained_deviance_all = [explained_deviance_test',...
    explained_deviance_tact_stim_onset_left_removed_test',...
    explained_deviance_tact_stim_onset_right_removed_test',...
    explained_deviance_delay_middle_left_removed_test',...
    explained_deviance_delay_middle_right_removed_test',...
    explained_deviance_left_lick_onset_removed_test',...
    explained_deviance_right_lick_onset_removed_test',...
    explained_deviance_reward_onset_removed_test',...
    explained_deviance_leftforepaw_movement_removed_test',...
    explained_deviance_rightforepaw_movement_removed_test',...
    explained_deviance_cell_coupling_removed_test'];

% Check if number of grouped task predictors is correct
number_task_predictor = size([size(resp_task1_var,2):size(explained_deviance_all,2)],2);
assert(number_task_predictor == 11)

% Note that explained_deviance idx of predictor groups is different from coeff idx
explained_deviance_all_reduction_temp = explained_deviance_all - explained_deviance_all(:,1); % For the task variable to be useful, explained deviance has to be smaller than the explained_deviance_test
explained_deviance_all_reduction = explained_deviance_all_reduction_temp(:,2:end); % Remove the full model.

%% Shuffle.
shuffle = 1000;
for shuffle_num = 1:1000
    disp(['Shuffle: ',num2str(shuffle_num)])
    clearvars -except shuffle combi multi_date all_session_num predictor_labels animal_ID date num_pca_dimensions fs_image Bpod_protocol shuffle task1_var resp_task1_var_test test_frame GLM_result shuffle_num explained_deviance_all explained_deviance_all_reduction_temp explained_deviance_all_reduction explained_deviance_all_shuffle explained_deviance_all_reduction_temp_shuffle explained_deviance_all_reduction_shuffle...
        tact_stim_onset_left_idx tact_stim_onset_right_idx delay_left_middle_idx delay_right_middle_idx left_lick_onset_idx right_lick_onset_idx reward_onset_idx leftforepaw_movement_idx rightforepaw_movement_idx
    
    % Full model.
    B0 = GLM_result.B0;
    coeff = GLM_result.coeff;
    predict_train = GLM_result.predict_train;
    predict_test = GLM_result.predict_test;
    y_test = GLM_result.y_test;
    y_hat_test = GLM_result.y_hat_test;
    y_null_test = GLM_result.y_null_test;
    L1_test = GLM_result.L1_test;
    L0_test = GLM_result.L0_test;
    LS_test = GLM_result.LS_test;
    explained_deviance_test = GLM_result.explained_deviance_test;
    
    % To look at individual task predictors in the delay period (Note that this is the reverse)
    delay_left_middle_1_idx = delay_left_middle_idx(1);
    delay_left_middle_2_idx = delay_left_middle_idx(2);
    delay_left_middle_3_idx = delay_left_middle_idx(3);
    delay_left_middle_4_idx = delay_left_middle_idx(4);
    delay_left_middle_5_idx = delay_left_middle_idx(5);
    
    delay_right_middle_1_idx = delay_right_middle_idx(1);
    delay_right_middle_2_idx = delay_right_middle_idx(2);
    delay_right_middle_3_idx = delay_right_middle_idx(3);
    delay_right_middle_4_idx = delay_right_middle_idx(4);
    delay_right_middle_5_idx = delay_right_middle_idx(5);
    
    % Mean y_hat variables
    mean_y_hat_tact_stim_onset_left = nan(length(GLM_result.B0),1);
    mean_y_hat_tact_stim_onset_right = nan(length(GLM_result.B0),1);
    mean_y_hat_delay_middle_left = nan(length(GLM_result.B0),1);
    mean_y_hat_delay_middle_right = nan(length(GLM_result.B0),1);
    % mean_y_hat_go_cue = nan(length(GLM_result.B0),1);
    mean_y_hat_left_lick_onset = nan(length(GLM_result.B0),1);
    mean_y_hat_right_lick_onset = nan(length(GLM_result.B0),1);
    mean_y_hat_reward_onset = nan(length(GLM_result.B0),1);
    mean_y_hat_leftforepaw_movement = nan(length(GLM_result.B0),1);
    mean_y_hat_rightforepaw_movement = nan(length(GLM_result.B0),1);
    
    mean_y_hat_delay_middle_left_1 = nan(length(GLM_result.B0),1);
    mean_y_hat_delay_middle_left_2 = nan(length(GLM_result.B0),1);
    mean_y_hat_delay_middle_left_3 = nan(length(GLM_result.B0),1);
    mean_y_hat_delay_middle_left_4 = nan(length(GLM_result.B0),1);
    mean_y_hat_delay_middle_right_1 = nan(length(GLM_result.B0),1);
    mean_y_hat_delay_middle_right_2 = nan(length(GLM_result.B0),1);
    mean_y_hat_delay_middle_right_3 = nan(length(GLM_result.B0),1);
    mean_y_hat_delay_middle_right_4 = nan(length(GLM_result.B0),1);
    
    % y_hat_variables
    y_hat_tact_stim_onset_left_removed = nan(length(GLM_result.B0),size(task1_var,1));
    y_hat_tact_stim_onset_left_removed_test= nan(length(GLM_result.B0),size(test_frame,2));
    L1_tact_stim_onset_left_removed_test = nan(length(GLM_result.B0),1);
    explained_deviance_tact_stim_onset_left_removed_test = nan(length(GLM_result.B0),1);
    
    y_hat_tact_stim_onset_right_removed = nan(length(GLM_result.B0),size(task1_var,1));
    y_hat_tact_stim_onset_right_removed_test= nan(length(GLM_result.B0),size(test_frame,2));
    L1_tact_stim_onset_right_removed_test = nan(length(GLM_result.B0),1);
    explained_deviance_tact_stim_onset_right_removed_test = nan(length(GLM_result.B0),1);
    
    y_hat_delay_middle_left_removed = nan(length(GLM_result.B0),size(task1_var,1));
    y_hat_delay_middle_left_removed_test= nan(length(GLM_result.B0),size(test_frame,2));
    L1_delay_middle_left_removed_test = nan(length(GLM_result.B0),1);
    explained_deviance_delay_middle_left_removed_test = nan(length(GLM_result.B0),1);
    
    y_hat_delay_middle_right_removed = nan(length(GLM_result.B0),size(task1_var,1));
    y_hat_delay_middle_right_removed_test= nan(length(GLM_result.B0),size(test_frame,2));
    L1_delay_middle_right_removed_test = nan(length(GLM_result.B0),1);
    explained_deviance_delay_middle_right_removed_test = nan(length(GLM_result.B0),1);
    
    y_hat_left_lick_onset_removed = nan(length(GLM_result.B0),size(task1_var,1));
    y_hat_left_lick_onset_removed_test= nan(length(GLM_result.B0),size(test_frame,2));
    L1_left_lick_onset_removed_test = nan(length(GLM_result.B0),1);
    explained_deviance_left_lick_onset_removed_test = nan(length(GLM_result.B0),1);
    
    y_hat_right_lick_onset_removed = nan(length(GLM_result.B0),size(task1_var,1));
    y_hat_right_lick_onset_removed_test= nan(length(GLM_result.B0),size(test_frame,2));
    L1_right_lick_onset_removed_test = nan(length(GLM_result.B0),1);
    explained_deviance_right_lick_onset_removed_test = nan(length(GLM_result.B0),1);
    
    y_hat_reward_onset_removed = nan(length(GLM_result.B0),size(task1_var,1));
    y_hat_reward_onset_removed_test= nan(length(GLM_result.B0),size(test_frame,2));
    L1_reward_onset_removed_test = nan(length(GLM_result.B0),1);
    explained_deviance_reward_onset_removed_test = nan(length(GLM_result.B0),1);
    
    y_hat_leftforepaw_movement_removed = nan(length(GLM_result.B0),size(task1_var,1));
    y_hat_leftforepaw_movement_removed_test= nan(length(GLM_result.B0),size(test_frame,2));
    L1_leftforepaw_movement_removed_test = nan(length(GLM_result.B0),1);
    explained_deviance_leftforepaw_movement_removed_test = nan(length(GLM_result.B0),1);
    
    y_hat_rightforepaw_movement_removed = nan(length(GLM_result.B0),size(task1_var,1));
    y_hat_rightforepaw_movement_removed_test = nan(length(GLM_result.B0),size(test_frame,2));
    L1_rightforepaw_movement_removed_test = nan(length(GLM_result.B0),1);
    explained_deviance_rightforepaw_movement_removed_test = nan(length(GLM_result.B0),1);
    
    y_hat_delay_middle_left_1_removed = nan(length(GLM_result.B0),size(task1_var,1));
    y_hat_delay_middle_left_1_removed_test = nan(length(GLM_result.B0),size(test_frame,2));
    L1_delay_middle_left_1_removed_test = nan(length(GLM_result.B0),1);
    explained_deviance_delay_middle_left_1_removed_test = nan(length(GLM_result.B0),1);
    
    y_hat_delay_middle_left_2_removed = nan(length(GLM_result.B0),size(task1_var,1));
    y_hat_delay_middle_left_2_removed_test = nan(length(GLM_result.B0),size(test_frame,2));
    L1_delay_middle_left_2_removed_test = nan(length(GLM_result.B0),1);
    explained_deviance_delay_middle_left_2_removed_test = nan(length(GLM_result.B0),1);
    
    y_hat_delay_middle_left_3_removed = nan(length(GLM_result.B0),size(task1_var,1));
    y_hat_delay_middle_left_3_removed_test = nan(length(GLM_result.B0),size(test_frame,2));
    L1_delay_middle_left_3_removed_test = nan(length(GLM_result.B0),1);
    explained_deviance_delay_middle_left_3_removed_test = nan(length(GLM_result.B0),1);
    
    y_hat_delay_middle_left_4_removed = nan(length(GLM_result.B0),size(task1_var,1));
    y_hat_delay_middle_left_4_removed_test = nan(length(GLM_result.B0),size(test_frame,2));
    L1_delay_middle_left_4_removed_test = nan(length(GLM_result.B0),1);
    explained_deviance_delay_middle_left_4_removed_test = nan(length(GLM_result.B0),1);
    
    y_hat_delay_middle_right_1_removed = nan(length(GLM_result.B0),size(task1_var,1));
    y_hat_delay_middle_right_1_removed_test = nan(length(GLM_result.B0),size(test_frame,2));
    L1_delay_middle_right_1_removed_test = nan(length(GLM_result.B0),1);
    explained_deviance_delay_middle_right_1_removed_test = nan(length(GLM_result.B0),1);
    
    y_hat_delay_middle_right_2_removed = nan(length(GLM_result.B0),size(task1_var,1));
    y_hat_delay_middle_right_2_removed_test = nan(length(GLM_result.B0),size(test_frame,2));
    L1_delay_middle_right_2_removed_test = nan(length(GLM_result.B0),1);
    explained_deviance_delay_middle_right_2_removed_test = nan(length(GLM_result.B0),1);
    
    y_hat_delay_middle_right_3_removed = nan(length(GLM_result.B0),size(task1_var,1));
    y_hat_delay_middle_right_3_removed_test = nan(length(GLM_result.B0),size(test_frame,2));
    L1_delay_middle_right_3_removed_test = nan(length(GLM_result.B0),1);
    explained_deviance_delay_middle_right_3_removed_test = nan(length(GLM_result.B0),1);
    
    y_hat_delay_middle_right_right_4_removed = nan(length(GLM_result.B0),size(task1_var,1));
    y_hat_delay_middle_right_4_removed_test = nan(length(GLM_result.B0),size(test_frame,2));
    L1_delay_middle_right_4_removed_test = nan(length(GLM_result.B0),1);
    explained_deviance_delay_middle_right_4_removed_test = nan(length(GLM_result.B0),1);
    
    % Initialise cell coupling variables using number of cells. Note that T1 and T2 coupling will be removed as a pair
    mean_y_hat_cell_coupling = nan(length(GLM_result.B0),length(GLM_result.B0)); % (Cell coupling predictors x predicted cells x 1)
    y_hat_cell_coupling_removed = nan(length(GLM_result.B0),length(GLM_result.B0),size(task1_var,1)); % (Cell coupling predictors x predicted cells x bins)
    y_hat_cell_coupling_removed_test = nan(length(GLM_result.B0),length(GLM_result.B0),size(test_frame,2));
    L1_cell_coupling_removed_test = nan(length(GLM_result.B0),length(GLM_result.B0),1);
    explained_deviance_cell_coupling_removed_test = nan(length(GLM_result.B0),length(GLM_result.B0),1);
    
    % Shuffle task1_var frames.
    % To be deterministic.
    rng(shuffle_num); % Seed.
    
    % Shuffle ~2 s long frame bins of task1_var.
    idx = 1:floor(size(task1_var,1)/round(2*fs_image))*round(2*fs_image); % Get a vector from 1 to the max product of round(2*fs_image).
    idx = reshape(idx,round(2*fs_image),[]);
    idx = reshape(idx(:,randperm(size(idx,2))),1,[]);
    idx = [idx,(length(idx) + 1):size(task1_var,1)];
    task1_var = task1_var(idx,:);
    
    % Get mean of cell coupling predictors using for loop
    for predictor_cell_num = 1 : max(unique(predictor_labels(~isnan(predictor_labels))))
        % Get predictor labels for t1 and t2
        predictor_idx = find(predictor_labels(~isnan(predictor_labels)) == predictor_cell_num) + max(find(isnan(predictor_labels)));
        % I.e. Taking the mean of the predicted activity at each predictor cell index
        mean_y_hat_cell_coupling(predictor_cell_num,:) = mean(exp(task1_var(:,predictor_idx)*coeff(:,predictor_idx)'));
    end
    
    % Scaling factor marginalizing out the effect of the other variables.
    % Getting the mean of the predicted activity
    mean_y_hat_tact_stim_onset_left = mean(exp(task1_var(:,tact_stim_onset_left_idx)*coeff(:,tact_stim_onset_left_idx)'));
    mean_y_hat_tact_stim_onset_right = mean(exp(task1_var(:,tact_stim_onset_right_idx)*coeff(:,tact_stim_onset_right_idx)'));
    mean_y_hat_delay_middle_left = mean(exp(task1_var(:,delay_left_middle_idx)*coeff(:,delay_left_middle_idx)'));
    mean_y_hat_delay_middle_right = mean(exp(task1_var(:,delay_right_middle_idx)*coeff(:,delay_right_middle_idx)'));
    
    % mean_y_hat_go_cue = mean(exp(task1_var(:,go_cue_idx)*coeff(:,go_cue_idx)'));
    mean_y_hat_left_lick_onset = mean(exp(task1_var(:,left_lick_onset_idx)*coeff(:,left_lick_onset_idx)'));
    mean_y_hat_right_lick_onset = mean(exp(task1_var(:,right_lick_onset_idx)*coeff(:,right_lick_onset_idx)'));
    mean_y_hat_reward_onset = mean(exp(task1_var(:,reward_onset_idx)*coeff(:,reward_onset_idx)'));
    mean_y_hat_leftforepaw_movement = mean(exp(task1_var(:,leftforepaw_movement_idx)*coeff(:,leftforepaw_movement_idx)'));
    mean_y_hat_rightforepaw_movement = mean(exp(task1_var(:,rightforepaw_movement_idx)*coeff(:,rightforepaw_movement_idx)'));
    
    mean_y_hat_delay_middle_left_1 = mean(exp(task1_var(:,delay_left_middle_1_idx)*coeff(:,delay_left_middle_1_idx)'));
    mean_y_hat_delay_middle_left_2 = mean(exp(task1_var(:,delay_left_middle_2_idx)*coeff(:,delay_left_middle_2_idx)'));
    mean_y_hat_delay_middle_left_3 = mean(exp(task1_var(:,delay_left_middle_3_idx)*coeff(:,delay_left_middle_3_idx)'));
    mean_y_hat_delay_middle_left_4 = mean(exp(task1_var(:,delay_left_middle_4_idx)*coeff(:,delay_left_middle_4_idx)'));
    mean_y_hat_delay_middle_right_1 = mean(exp(task1_var(:,delay_right_middle_1_idx)*coeff(:,delay_right_middle_1_idx)'));
    mean_y_hat_delay_middle_right_2 = mean(exp(task1_var(:,delay_right_middle_2_idx)*coeff(:,delay_right_middle_2_idx)'));
    mean_y_hat_delay_middle_right_3 = mean(exp(task1_var(:,delay_right_middle_3_idx)*coeff(:,delay_right_middle_3_idx)'));
    mean_y_hat_delay_middle_right_4 = mean(exp(task1_var(:,delay_right_middle_4_idx)*coeff(:,delay_right_middle_4_idx)'));
    
    clear y_test y_nul_test L0_test LS_test
    for cell_num = 1: length(GLM_result.B0)
        % Model assessment by measuring explained deviance based on Benjamin et al 2018.
        y_test(cell_num,:) = resp_task1_var_test(:,cell_num);
        y_null_test(cell_num) = mean(resp_task1_var_test(:,cell_num),1);
        L0_test(cell_num) = sum(y_test(cell_num,:).*log(eps + y_null_test(cell_num)) - y_null_test(cell_num)); % Null model.
        LS_test(cell_num) = sum(y_test(cell_num,:).*log(eps + y_test(cell_num,:)) - y_test(cell_num,:)); % Saturated model.
        
    end
    
    % Obtain partial models.
    %% Get predicted activity after removing the specific predictor
    
    assert( max(find(isnan(predictor_labels))) == max(rightforepaw_movement_idx));
    start_coupling_idx = max(find(isnan(predictor_labels))) + 1;
    end_coupling_idx = size(task1_var,2);
    
    % No variable removed
    matrix_y_hat_no_var_removed = ...
        exp(task1_var(:,start_coupling_idx:end_coupling_idx)*coeff(:,start_coupling_idx:end_coupling_idx)').*...
        exp(task1_var(:,tact_stim_onset_left_idx)*coeff(:,tact_stim_onset_left_idx)').*...
        exp(task1_var(:,tact_stim_onset_right_idx)*coeff(:,tact_stim_onset_right_idx)').*...
        exp(task1_var(:,delay_left_middle_idx)*coeff(:,delay_left_middle_idx)').*...
        exp(task1_var(:,delay_right_middle_idx)*coeff(:,delay_right_middle_idx)').*...
        exp(task1_var(:,left_lick_onset_idx)*coeff(:,left_lick_onset_idx)').*...
        exp(task1_var(:,right_lick_onset_idx)*coeff(:,right_lick_onset_idx)').*...
        exp(task1_var(:,reward_onset_idx)*coeff(:,reward_onset_idx)').*...
        exp(task1_var(:,leftforepaw_movement_idx)*coeff(:,leftforepaw_movement_idx)').*...
        exp(task1_var(:,rightforepaw_movement_idx)*coeff(:,rightforepaw_movement_idx)').*...
        exp(B0);
    y_hat_no_var_removed = matrix_y_hat_no_var_removed';
    
    % Remove cell_coupling using (Predictor cells (Activity removed) x Predicted cells x Bins) size for y_hat_cell_coupling_removed
    parfor predictor_cell_num = 1 : max(unique(predictor_labels(~isnan(predictor_labels))))
        predictor_idx = find(predictor_labels(~isnan(predictor_labels)) == predictor_cell_num) + max(find(isnan(predictor_labels)));
        % Get the index of all cell coupling predictors except the one being removed
        all_predictors_except = setdiff(start_coupling_idx:end_coupling_idx, predictor_idx)
        
        % Get 3D tensor of y_hat, with the first dim. being the predictor cell that is removed
        y_hat_cell_coupling_removed(predictor_cell_num,:,:) = ...
            (exp(task1_var(:,all_predictors_except)*coeff(:,all_predictors_except)').*...
            mean_y_hat_cell_coupling(predictor_cell_num,:).*...
            exp(task1_var(:,tact_stim_onset_left_idx)*coeff(:,tact_stim_onset_left_idx)').*...
            exp(task1_var(:,tact_stim_onset_right_idx)*coeff(:,tact_stim_onset_right_idx)').*...
            exp(task1_var(:,delay_left_middle_idx)*coeff(:,delay_left_middle_idx)').*...
            exp(task1_var(:,delay_right_middle_idx)*coeff(:,delay_right_middle_idx)').*...
            exp(task1_var(:,left_lick_onset_idx)*coeff(:,left_lick_onset_idx)').*...
            exp(task1_var(:,right_lick_onset_idx)*coeff(:,right_lick_onset_idx)').*...
            exp(task1_var(:,reward_onset_idx)*coeff(:,reward_onset_idx)').*...
            exp(task1_var(:,leftforepaw_movement_idx)*coeff(:,leftforepaw_movement_idx)').*...
            exp(task1_var(:,rightforepaw_movement_idx)*coeff(:,rightforepaw_movement_idx)').*...
            exp(B0))';
    end
    
    % Remove tactile stimulus onset left.
    matrix_y_hat_tact_stim_onset_left_removed = ...
        mean_y_hat_tact_stim_onset_left.*...
        exp(task1_var(:,start_coupling_idx:end_coupling_idx)*coeff(:,start_coupling_idx:end_coupling_idx)').*...
        exp(task1_var(:,tact_stim_onset_right_idx)*coeff(:,tact_stim_onset_right_idx)').*...
        exp(task1_var(:,delay_left_middle_idx)*coeff(:,delay_left_middle_idx)').*...
        exp(task1_var(:,delay_right_middle_idx)*coeff(:,delay_right_middle_idx)').*...
        exp(task1_var(:,left_lick_onset_idx)*coeff(:,left_lick_onset_idx)').*...
        exp(task1_var(:,right_lick_onset_idx)*coeff(:,right_lick_onset_idx)').*...
        exp(task1_var(:,reward_onset_idx)*coeff(:,reward_onset_idx)').*...
        exp(task1_var(:,leftforepaw_movement_idx)*coeff(:,leftforepaw_movement_idx)').*...
        exp(task1_var(:,rightforepaw_movement_idx)*coeff(:,rightforepaw_movement_idx)').*...
        exp(B0);
    y_hat_tact_stim_onset_left_removed = matrix_y_hat_tact_stim_onset_left_removed';
    
    % Remove tactile stimulus onset right.
    matrix_y_hat_tact_stim_onset_right_removed = ...
        exp(task1_var(:,start_coupling_idx:end_coupling_idx)*coeff(:,start_coupling_idx:end_coupling_idx)').*...
        exp(task1_var(:,tact_stim_onset_left_idx)*coeff(:,tact_stim_onset_left_idx)').*...
        mean_y_hat_tact_stim_onset_right.*...
        exp(task1_var(:,delay_left_middle_idx)*coeff(:,delay_left_middle_idx)').*...
        exp(task1_var(:,delay_right_middle_idx)*coeff(:,delay_right_middle_idx)').*...
        exp(task1_var(:,left_lick_onset_idx)*coeff(:,left_lick_onset_idx)').*...
        exp(task1_var(:,right_lick_onset_idx)*coeff(:,right_lick_onset_idx)').*...
        exp(task1_var(:,reward_onset_idx)*coeff(:,reward_onset_idx)').*...
        exp(task1_var(:,leftforepaw_movement_idx)*coeff(:,leftforepaw_movement_idx)').*...
        exp(task1_var(:,rightforepaw_movement_idx)*coeff(:,rightforepaw_movement_idx)').*...
        exp(B0);
    y_hat_tact_stim_onset_right_removed = matrix_y_hat_tact_stim_onset_right_removed';
    
    % Remove delay middle left.
    matrix_y_hat_delay_middle_left_removed = ...
        exp(task1_var(:,start_coupling_idx:end_coupling_idx)*coeff(:,start_coupling_idx:end_coupling_idx)').*...
        exp(task1_var(:,tact_stim_onset_left_idx)*coeff(:,tact_stim_onset_left_idx)').*...
        exp(task1_var(:,tact_stim_onset_right_idx)*coeff(:,tact_stim_onset_right_idx)').*...
        mean_y_hat_delay_middle_left.*...
        exp(task1_var(:,delay_right_middle_idx)*coeff(:,delay_right_middle_idx)').*...
        exp(task1_var(:,left_lick_onset_idx)*coeff(:,left_lick_onset_idx)').*...
        exp(task1_var(:,right_lick_onset_idx)*coeff(:,right_lick_onset_idx)').*...
        exp(task1_var(:,reward_onset_idx)*coeff(:,reward_onset_idx)').*...
        exp(task1_var(:,leftforepaw_movement_idx)*coeff(:,leftforepaw_movement_idx)').*...
        exp(task1_var(:,rightforepaw_movement_idx)*coeff(:,rightforepaw_movement_idx)').*...
        exp(B0);
    y_hat_delay_middle_left_removed = matrix_y_hat_delay_middle_left_removed';
    
    % Remove delay middle right
    matrix_y_hat_delay_middle_right_removed = ...
        exp(task1_var(:,start_coupling_idx:end_coupling_idx)*coeff(:,start_coupling_idx:end_coupling_idx)').*...
        exp(task1_var(:,tact_stim_onset_left_idx)*coeff(:,tact_stim_onset_left_idx)').*...
        exp(task1_var(:,tact_stim_onset_right_idx)*coeff(:,tact_stim_onset_right_idx)').*...
        exp(task1_var(:,delay_left_middle_idx)*coeff(:,delay_left_middle_idx)').*...
        mean_y_hat_delay_middle_right.*...
        exp(task1_var(:,left_lick_onset_idx)*coeff(:,left_lick_onset_idx)').*...
        exp(task1_var(:,right_lick_onset_idx)*coeff(:,right_lick_onset_idx)').*...
        exp(task1_var(:,reward_onset_idx)*coeff(:,reward_onset_idx)').*...
        exp(task1_var(:,leftforepaw_movement_idx)*coeff(:,leftforepaw_movement_idx)').*...
        exp(task1_var(:,rightforepaw_movement_idx)*coeff(:,rightforepaw_movement_idx)').*...
        exp(B0);
    y_hat_delay_middle_right_removed = matrix_y_hat_delay_middle_right_removed';
    
    % Remove left lick onset.
    matrix_y_hat_left_lick_onset_removed = ...
        exp(task1_var(:,start_coupling_idx:end_coupling_idx)*coeff(:,start_coupling_idx:end_coupling_idx)').*...
        exp(task1_var(:,tact_stim_onset_left_idx)*coeff(:,tact_stim_onset_left_idx)').*...
        exp(task1_var(:,tact_stim_onset_right_idx)*coeff(:,tact_stim_onset_right_idx)').*...
        exp(task1_var(:,delay_left_middle_idx)*coeff(:,delay_left_middle_idx)').*...
        exp(task1_var(:,delay_right_middle_idx)*coeff(:,delay_right_middle_idx)').*...
        mean_y_hat_left_lick_onset.*...
        exp(task1_var(:,right_lick_onset_idx)*coeff(:,right_lick_onset_idx)').*...
        exp(task1_var(:,reward_onset_idx)*coeff(:,reward_onset_idx)').*...
        exp(task1_var(:,leftforepaw_movement_idx)*coeff(:,leftforepaw_movement_idx)').*...
        exp(task1_var(:,rightforepaw_movement_idx)*coeff(:,rightforepaw_movement_idx)').*...
        exp(B0);
    y_hat_left_lick_onset_removed = matrix_y_hat_left_lick_onset_removed';
    
    % Remove right lick onset.
    matrix_y_hat_right_lick_onset_removed = ...
        exp(task1_var(:,start_coupling_idx:end_coupling_idx)*coeff(:,start_coupling_idx:end_coupling_idx)').*...
        exp(task1_var(:,tact_stim_onset_left_idx)*coeff(:,tact_stim_onset_left_idx)').*...
        exp(task1_var(:,tact_stim_onset_right_idx)*coeff(:,tact_stim_onset_right_idx)').*...
        exp(task1_var(:,delay_left_middle_idx)*coeff(:,delay_left_middle_idx)').*...
        exp(task1_var(:,delay_right_middle_idx)*coeff(:,delay_right_middle_idx)').*...
        exp(task1_var(:,left_lick_onset_idx)*coeff(:,left_lick_onset_idx)').*...
        mean_y_hat_right_lick_onset.*...
        exp(task1_var(:,reward_onset_idx)*coeff(:,reward_onset_idx)').*...
        exp(task1_var(:,leftforepaw_movement_idx)*coeff(:,leftforepaw_movement_idx)').*...
        exp(task1_var(:,rightforepaw_movement_idx)*coeff(:,rightforepaw_movement_idx)').*...
        exp(B0);
    y_hat_right_lick_onset_removed = matrix_y_hat_right_lick_onset_removed';
    
    % Remove reward onset.
    matrix_y_hat_reward_onset_removed = ...
        exp(task1_var(:,start_coupling_idx:end_coupling_idx)*coeff(:,start_coupling_idx:end_coupling_idx)').*...
        exp(task1_var(:,tact_stim_onset_left_idx)*coeff(:,tact_stim_onset_left_idx)').*...
        exp(task1_var(:,tact_stim_onset_right_idx)*coeff(:,tact_stim_onset_right_idx)').*...
        exp(task1_var(:,delay_left_middle_idx)*coeff(:,delay_left_middle_idx)').*...
        exp(task1_var(:,delay_right_middle_idx)*coeff(:,delay_right_middle_idx)').*...
        exp(task1_var(:,left_lick_onset_idx)*coeff(:,left_lick_onset_idx)').*...
        exp(task1_var(:,right_lick_onset_idx)*coeff(:,right_lick_onset_idx)').*...
        mean_y_hat_reward_onset.*...
        exp(task1_var(:,leftforepaw_movement_idx)*coeff(:,leftforepaw_movement_idx)').*...
        exp(task1_var(:,rightforepaw_movement_idx)*coeff(:,rightforepaw_movement_idx)').*...
        exp(B0);
    y_hat_reward_onset_removed = matrix_y_hat_reward_onset_removed';
    
    % Remove left forepaw movement.
    matrix_y_hat_leftforepaw_movement_removed = ...
        exp(task1_var(:,start_coupling_idx:end_coupling_idx)*coeff(:,start_coupling_idx:end_coupling_idx)').*...
        exp(task1_var(:,tact_stim_onset_left_idx)*coeff(:,tact_stim_onset_left_idx)').*...
        exp(task1_var(:,tact_stim_onset_right_idx)*coeff(:,tact_stim_onset_right_idx)').*...
        exp(task1_var(:,delay_left_middle_idx)*coeff(:,delay_left_middle_idx)').*...
        exp(task1_var(:,delay_right_middle_idx)*coeff(:,delay_right_middle_idx)').*...
        exp(task1_var(:,left_lick_onset_idx)*coeff(:,left_lick_onset_idx)').*...
        exp(task1_var(:,right_lick_onset_idx)*coeff(:,right_lick_onset_idx)').*...
        exp(task1_var(:,reward_onset_idx)*coeff(:,reward_onset_idx)').*...
        mean_y_hat_leftforepaw_movement.*...
        exp(task1_var(:,rightforepaw_movement_idx)*coeff(:,rightforepaw_movement_idx)').*...
        exp(B0);
    y_hat_leftforepaw_movement_removed  = matrix_y_hat_leftforepaw_movement_removed';
    
    % Remove right forepaw movement.
    matrix_y_hat_rightforepaw_movement_removed = ...
        exp(task1_var(:,start_coupling_idx:end_coupling_idx)*coeff(:,start_coupling_idx:end_coupling_idx)').*...
        exp(task1_var(:,tact_stim_onset_left_idx)*coeff(:,tact_stim_onset_left_idx)').*...
        exp(task1_var(:,tact_stim_onset_right_idx)*coeff(:,tact_stim_onset_right_idx)').*...
        exp(task1_var(:,delay_left_middle_idx)*coeff(:,delay_left_middle_idx)').*...
        exp(task1_var(:,delay_right_middle_idx)*coeff(:,delay_right_middle_idx)').*...
        exp(task1_var(:,left_lick_onset_idx)*coeff(:,left_lick_onset_idx)').*...
        exp(task1_var(:,right_lick_onset_idx)*coeff(:,right_lick_onset_idx)').*...
        exp(task1_var(:,reward_onset_idx)*coeff(:,reward_onset_idx)').*...
        exp(task1_var(:,leftforepaw_movement_idx)*coeff(:,leftforepaw_movement_idx)').*...
        mean_y_hat_rightforepaw_movement.*...
        exp(B0);
    y_hat_rightforepaw_movement_removed  = matrix_y_hat_rightforepaw_movement_removed';
    
    %% Get explained deviance of removed model
    y_hat_no_var_removed_test = y_hat_no_var_removed(:,test_frame);
    L1_no_var_removed_test = sum(y_test.*log(eps + y_hat_no_var_removed_test) - y_hat_no_var_removed_test,2); % Partial model.
    explained_deviance_no_var_removed_test = 1 - (LS_test - L1_no_var_removed_test')./(LS_test - L0_test);
    
    clear explained_deviance_cell_coupling_removed_test
    parfor predictor_cell_num = 1 : max(unique(predictor_labels(~isnan(predictor_labels))))
        y_hat_cell_coupling_removed_test = squeeze(y_hat_cell_coupling_removed(predictor_cell_num, :, test_frame)); % Only take into account the tests trials
        L1_cell_coupling_removed_test = sum(y_test.*log(eps + y_hat_cell_coupling_removed_test) - y_hat_cell_coupling_removed_test,2); % Partial model.
        
        % Arrange predictor cell_num removed idx in columns
        explained_deviance_cell_coupling_removed_test(predictor_cell_num,:) = 1 - (LS_test - L1_cell_coupling_removed_test')./(LS_test - L0_test);
    end
    
    y_hat_tact_stim_onset_left_removed_test = y_hat_tact_stim_onset_left_removed(:,test_frame); % Only take into account the tests trials
    L1_tact_stim_onset_left_removed_test = sum(y_test.*log(eps + y_hat_tact_stim_onset_left_removed_test) - y_hat_tact_stim_onset_left_removed_test,2); % Partial model.
    explained_deviance_tact_stim_onset_left_removed_test = 1 - (LS_test - L1_tact_stim_onset_left_removed_test')./(LS_test - L0_test);
    
    y_hat_tact_stim_onset_right_removed_test = y_hat_tact_stim_onset_right_removed(:,test_frame);
    L1_tact_stim_onset_right_removed_test = sum(y_test.*log(eps + y_hat_tact_stim_onset_right_removed_test) - y_hat_tact_stim_onset_right_removed_test,2); % Partial model.
    explained_deviance_tact_stim_onset_right_removed_test = 1 - (LS_test - L1_tact_stim_onset_right_removed_test')./(LS_test - L0_test);
    
    y_hat_delay_middle_left_removed_test = y_hat_delay_middle_left_removed(:,test_frame);
    L1_delay_middle_left_removed_test = sum(y_test.*log(eps + y_hat_delay_middle_left_removed_test) - y_hat_delay_middle_left_removed_test,2); % Partial model.
    explained_deviance_delay_middle_left_removed_test = 1 - (LS_test - L1_delay_middle_left_removed_test')./(LS_test - L0_test);
    
    y_hat_delay_middle_right_removed_test = y_hat_delay_middle_right_removed(:,test_frame);
    L1_delay_middle_right_removed_test = sum(y_test.*log(eps + y_hat_delay_middle_right_removed_test) - y_hat_delay_middle_right_removed_test,2); % Partial model.
    explained_deviance_delay_middle_right_removed_test = 1 - (LS_test - L1_delay_middle_right_removed_test')./(LS_test - L0_test);
    
    y_hat_left_lick_onset_removed_test = y_hat_left_lick_onset_removed(:,test_frame);
    L1_left_lick_onset_removed_test = sum(y_test.*log(eps + y_hat_left_lick_onset_removed_test) - y_hat_left_lick_onset_removed_test,2); % Partial model.
    explained_deviance_left_lick_onset_removed_test = 1 - (LS_test - L1_left_lick_onset_removed_test')./(LS_test - L0_test);
    
    y_hat_right_lick_onset_removed_test = y_hat_right_lick_onset_removed(:,test_frame);
    L1_right_lick_onset_removed_test = sum(y_test.*log(eps + y_hat_right_lick_onset_removed_test) - y_hat_right_lick_onset_removed_test,2); % Partial model.
    explained_deviance_right_lick_onset_removed_test = 1 - (LS_test - L1_right_lick_onset_removed_test')./(LS_test - L0_test);
    
    y_hat_reward_onset_removed_test = y_hat_reward_onset_removed(:,test_frame);
    L1_reward_onset_removed_test = sum(y_test.*log(eps + y_hat_reward_onset_removed_test) - y_hat_reward_onset_removed_test,2); % Partial model.
    explained_deviance_reward_onset_removed_test = 1 - (LS_test - L1_reward_onset_removed_test')./(LS_test - L0_test);
    
    y_hat_leftforepaw_movement_removed_test = y_hat_leftforepaw_movement_removed(:,test_frame);
    L1_leftforepaw_movement_removed_test = sum(y_test.*log(eps + y_hat_leftforepaw_movement_removed_test) - y_hat_leftforepaw_movement_removed_test,2); % Partial model.
    explained_deviance_leftforepaw_movement_removed_test = 1 - (LS_test - L1_leftforepaw_movement_removed_test')./(LS_test - L0_test);
    
    y_hat_rightforepaw_movement_removed_test = y_hat_rightforepaw_movement_removed(:,test_frame);
    L1_rightforepaw_movement_removed_test = sum(y_test.*log(eps + y_hat_rightforepaw_movement_removed_test) - y_hat_rightforepaw_movement_removed_test,2); % Partial model.
    explained_deviance_rightforepaw_movement_removed_test = 1 - (LS_test - L1_rightforepaw_movement_removed_test')./(LS_test - L0_test);
    
    % Assess the partial models. (Cells x Number of removed predictor groups)
    explained_deviance_all_shuffle{shuffle_num} = [explained_deviance_no_var_removed_test',...
        explained_deviance_tact_stim_onset_left_removed_test',...
        explained_deviance_tact_stim_onset_right_removed_test',...
        explained_deviance_delay_middle_left_removed_test',...
        explained_deviance_delay_middle_right_removed_test',...
        explained_deviance_left_lick_onset_removed_test',...
        explained_deviance_right_lick_onset_removed_test',...
        explained_deviance_reward_onset_removed_test',...
        explained_deviance_leftforepaw_movement_removed_test',...
        explained_deviance_rightforepaw_movement_removed_test',...
        explained_deviance_cell_coupling_removed_test'];
    
    explained_deviance_all_reduction_temp_shuffle{shuffle_num} = explained_deviance_all_shuffle{shuffle_num} - explained_deviance_all_shuffle{shuffle_num}(:,1);
    explained_deviance_all_reduction_shuffle{shuffle_num} = explained_deviance_all_reduction_temp_shuffle{shuffle_num}(:,2:end); % Remove the full model.
end

%% Concatenate.
explained_deviance_all_shuffle_concat = [];
explained_deviance_all_reduction_shuffle_concat = [];
for shuffle_num = 1:shuffle
    explained_deviance_all_shuffle_concat = cat(3,explained_deviance_all_shuffle_concat,explained_deviance_all_shuffle{shuffle_num});
    explained_deviance_all_reduction_shuffle_concat = cat(3,explained_deviance_all_reduction_shuffle_concat,explained_deviance_all_reduction_shuffle{shuffle_num});
end

% Obtain p-values.
for cell_num = 1:length(GLM_result.B0)
    % Compare actual explained deviance against shuffled explained deviance. Out of 1000 shuffles, how many points > 0? The number of points > 0 can be used to calculate the p value
    p_value(cell_num) = sum((squeeze(explained_deviance_all_shuffle_concat(cell_num,1,:)) - explained_deviance_all(cell_num,1)) > 0)/shuffle_num;
    
    % Compare substracted explained deviance against shuffled substracted explained deviance. Out of 1000 shuffles, how many points < 0? The number of points < 0 can be used to calculate the p value
    p_value_task1_var(cell_num,:) = sum((squeeze(explained_deviance_all_reduction_shuffle_concat(cell_num,:,:))' - explained_deviance_all_reduction(cell_num,:)) < 0)/shuffle_num;
end

% Save the contribution matrix.
save('GLM_significance.mat','explained_deviance_all_reduction','p_value','p_value_task1_var')