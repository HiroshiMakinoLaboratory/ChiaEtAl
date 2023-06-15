%% Get predicted activity from expert session
function ablated_data = get_ablated_data(parameters,data)
fs_image = parameters.fs_image;
fixed_num_shuffle = 100; % Number of shuffles
epoch_start = round(fs_image*3); % +3s from start of (-4s) ITI
epoch_end = round(fs_image*8); % 1s ITI, 1s stim, 2s Delay, 1s action
region_idx = [1 6 7 8]; % ALM, vS1, RSC, PPC

% Load dataset containing raw activity
all_norm_tact_stim_onset_left_convolved_trial = data.predictors.all_norm_tact_stim_onset_left_convolved_trial;
all_norm_tact_stim_onset_right_convolved_trial = data.predictors.all_norm_tact_stim_onset_right_convolved_trial;
all_norm_delay_middle_left_convolved_trial = data.predictors.all_norm_delay_middle_left_convolved_trial;
all_norm_delay_middle_right_convolved_trial = data.predictors.all_norm_delay_middle_right_convolved_trial;
all_norm_lick_left_onset_convolved_trial = data.predictors.all_norm_lick_left_onset_convolved_trial;
all_norm_lick_right_onset_convolved_trial = data.predictors.all_norm_lick_right_onset_convolved_trial;
all_norm_reward_onset_convolved_trial = data.predictors.all_norm_reward_onset_convolved_trial;
all_norm_leftforepaw_movement_convolved_trial = data.predictors.all_norm_leftforepaw_movement_convolved_trial;
all_norm_rightforepaw_movement_convolved_trial = data.predictors.all_norm_rightforepaw_movement_convolved_trial;
all_norm_resp_task1_var_t1_trial = data.coupling_activity.all_norm_resp_task1_var_t1_trial;
all_norm_resp_task1_var_t2_trial = data.coupling_activity.all_norm_resp_task1_var_t2_trial;
all_resp_var_trial_concat = data.activity.all_resp_var_trial_concat;
all_trial_type = data.trial_idx.all_trial_type;
all_valid_cells = data.cell_idx.all_valid_cells;
all_region_labels = data.cell_idx.all_region_labels;
all_predictor_labels = data.cell_idx.all_predictor_labels;
coupling_sig_cells = data.cell_idx.coupling_sig_cells;
all_coeff = data.cell_idx.all_coeff;
all_GLM = data.all_GLM;

%% Part 1: Reconstruct activity after ablating specific functional coupling
for animal_num = 1:14    
    task_variable_sig_cell{1} =  all_valid_cells;
    for session_num = 1
        % Concatenate all trials for task variables
        clear task1_var_temp left_trial_idx right_trial_idx
        % Get Task1_var_temp which are the normalised predictors (i.e.task + coupling predictors)
        % Concat only delay epoch together
        for trial_num = 1:size( all_norm_tact_stim_onset_left_convolved_trial{animal_num}{session_num},2)
            task1_var_temp{animal_num}{session_num}{trial_num} = [
                all_norm_tact_stim_onset_left_convolved_trial{animal_num}{session_num}{trial_num}(epoch_start:epoch_end,:),...
                all_norm_tact_stim_onset_right_convolved_trial{animal_num}{session_num}{trial_num}(epoch_start:epoch_end,:),...
                all_norm_delay_middle_left_convolved_trial{animal_num}{session_num}{trial_num}(epoch_start:epoch_end,:),...
                all_norm_delay_middle_right_convolved_trial{animal_num}{session_num}{trial_num}(epoch_start:epoch_end,:),...
                all_norm_lick_left_onset_convolved_trial{animal_num}{session_num}{trial_num}(epoch_start:epoch_end,:),...
                all_norm_lick_right_onset_convolved_trial{animal_num}{session_num}{trial_num}(epoch_start:epoch_end,:),...
                all_norm_reward_onset_convolved_trial{animal_num}{session_num}{trial_num}(epoch_start:epoch_end,:),...
                all_norm_leftforepaw_movement_convolved_trial{animal_num}{session_num}{trial_num}(epoch_start:epoch_end,:),...
                all_norm_rightforepaw_movement_convolved_trial{animal_num}{session_num}{trial_num}(epoch_start:epoch_end,:)];
            
            % Reorganise the t1 and t2 predictors and concatenate them into task1_var_temp
            for predictor_num = 1: size(all_norm_resp_task1_var_t1_trial{animal_num}{session_num}{trial_num},2)
                task1_var_temp{animal_num}{session_num}{trial_num} = [task1_var_temp{animal_num}{session_num}{trial_num},all_norm_resp_task1_var_t1_trial{animal_num}{session_num}{trial_num}(epoch_start:epoch_end,predictor_num)];
                task1_var_temp{animal_num}{session_num}{trial_num} = [task1_var_temp{animal_num}{session_num}{trial_num},all_norm_resp_task1_var_t2_trial{animal_num}{session_num}{trial_num}(epoch_start:epoch_end,predictor_num)];
            end
        end
        
        %% Stratifying trial types (correct left, incorrect left, etc)
        % Get trial types
        left_trial_idx = sort([find(all_trial_type{animal_num}{session_num} == 1)])';
        right_trial_idx = sort([find(all_trial_type{animal_num}{session_num} == 3)])';
        
        clear left_all_resp_var_trial_concat left_all_task1_var right_all_resp_var_trial_concat right_all_task1_var
        
        % Right choices
        right_all_resp_var_trial_concat{animal_num}{session_num} = [];
        right_all_task1_var{animal_num}{session_num} = [];
        for trial_num = right_trial_idx
            right_all_resp_var_trial_concat{animal_num}{session_num} = [right_all_resp_var_trial_concat{animal_num}{session_num} all_resp_var_trial_concat{animal_num}{session_num}{trial_num}(epoch_start:epoch_end,:)'];
            right_all_task1_var{animal_num}{session_num}  = [right_all_task1_var{animal_num}{session_num} task1_var_temp{animal_num}{session_num}{trial_num}'];
        end
        
        % left choices
        left_all_resp_var_trial_concat{animal_num}{session_num} = [];
        left_all_task1_var{animal_num}{session_num} = [];
        for trial_num = left_trial_idx
            left_all_resp_var_trial_concat{animal_num}{session_num} = [left_all_resp_var_trial_concat{animal_num}{session_num} all_resp_var_trial_concat{animal_num}{session_num}{trial_num}(epoch_start:epoch_end,:)'];
            left_all_task1_var{animal_num}{session_num}  = [left_all_task1_var{animal_num}{session_num} task1_var_temp{animal_num}{session_num}{trial_num}'];
        end
        
        %% Get predicted activity using coeff
        B0 = all_GLM{animal_num}{session_num}.B0;
        for region_num = region_idx
            other_region_idx = region_idx;
            for other_region_num = other_region_idx
                for i = 1:length(intersect(task_variable_sig_cell{1}{animal_num}{session_num},find(all_region_labels{animal_num}{session_num}==region_num))')
                    cell_idx = intersect(task_variable_sig_cell{1}{animal_num}{session_num},find(all_region_labels{animal_num}{session_num}==region_num))';
                    cell_num = cell_idx(i);
                    
                    % Intialise the variables
                    control_cell_right_predictors_predicted_activity{animal_num}{session_num}{region_num,other_region_num}{i} = [];
                    original_cell_right_predictors_predicted_activity{animal_num}{session_num}{region_num,other_region_num}{i} = [];
                    control_cell_left_predictors_predicted_activity{animal_num}{session_num}{region_num,other_region_num}{i} = [];
                    original_cell_left_predictors_predicted_activity{animal_num}{session_num}{region_num,other_region_num}{i} = [];
                end
            end
        end
        
        % Get predicted activity by reconstructing the GLM
        % Reconstruct control vs silenced population
        for region_num = region_idx
            for other_region_num = region_idx
                for i = 1:length(intersect(task_variable_sig_cell{1}{animal_num}{session_num},find(all_region_labels{animal_num}{session_num}==region_num))')
                    cell_idx = intersect(task_variable_sig_cell{1}{animal_num}{session_num},find(all_region_labels{animal_num}{session_num}==region_num))';
                    cell_num = cell_idx(i);
                    
                    % Get animal specific predictor labels
                    predictor_labels = all_predictor_labels{animal_num}{session_num};
                    predictor_idx = find(~isnan(predictor_labels));
                    
                    % Adjust the predictor idx to remove the same number of cells across diff regions
                    region_predictor_idx = find(ismember(predictor_labels,(find(all_region_labels{animal_num}{session_num} == other_region_num))));
                    
                    % Selectively remove fixed number of predictors
                    region_cells = intersect(find(all_region_labels{animal_num}{session_num} == other_region_num), task_variable_sig_cell{1}{animal_num}{session_num} );
                    non_region_cells = intersect(find(all_region_labels{animal_num}{session_num} ~= other_region_num), task_variable_sig_cell{1}{animal_num}{session_num} );
                    
                    clear   region_task_sig_cell to_remove_cells randomly_shuffled_to_remove_region_predictor_idx removed_predictor_idx control_to_remove_predictor_idx non_region_task_cells...
                        control_to_remove_cells control_shuffled_to_remove_region_predictor_idx control_predictor_idx non_region_task_cells
                    removed_predictor_idx = [];
                    control_predictor_idx = [];
                    
                    % These are cells that are 1. other region specific and 2. significantly coupled
                    region_task_sig_cell =intersect(region_cells, coupling_sig_cells{animal_num}{session_num}{cell_num} )';
                    non_region_task_sig_cell =intersect(non_region_cells, coupling_sig_cells{animal_num}{session_num}{cell_num} )';
                    all_region_task_sig_cell = intersect(task_variable_sig_cell{1}{animal_num}{session_num}, coupling_sig_cells{animal_num}{session_num}{cell_num} );
                    if ~isempty(region_task_sig_cell)
                        clear size_combi_vector num_cell_removal
                        
                        for hide = 1
                            num_cell_removal = round(length(region_task_sig_cell));
                            rng(2020);
                            % Each shuffle is to subsample the same number of cells from outside region or within region
                            for shuffle_num_cell = 1:fixed_num_shuffle
                                % Fix same number of cells
                                min_num_cells = min([length(region_task_sig_cell) length(non_region_task_sig_cell) ]);
                                rand_selected_region_task_sig_cell = region_task_sig_cell(randsample(1:length(region_task_sig_cell),min_num_cells));
                                
                                % Randomly select using the min number of cells
                                rand_selected_all_region_task_sig_cell = all_region_task_sig_cell(randsample(1:length(all_region_task_sig_cell),min_num_cells));
                                rand_selected_non_region_task_sig_cell = non_region_task_sig_cell(randsample(1:length(non_region_task_sig_cell),min_num_cells));
                                
                                % All the sig coupling n region cells predictor (CONTROL)
                                sig_coupling_predictor_idx = find(ismember(predictor_labels ,rand_selected_region_task_sig_cell));
                                
                                % All the sig coupling n non-region cells (ABLATED)
                                sig_outside_coupling_predictor_idx = find(ismember(predictor_labels ,rand_selected_non_region_task_sig_cell));
                                
                                % Reconstruct right trials
                                if ~isempty(region_cells) & ~isempty(right_all_task1_var{animal_num}{session_num})
                                    control_cell_right_predictors_predicted_activity{animal_num}{session_num}{region_num,other_region_num}{i}(:,shuffle_num_cell) =  exp(right_all_task1_var{animal_num}{session_num}(sig_outside_coupling_predictor_idx,:)' * all_coeff{animal_num}{session_num}(cell_num,sig_outside_coupling_predictor_idx)')... % Y hat for useful cell predictors
                                        .* exp(B0(cell_num));
                                    original_cell_right_predictors_predicted_activity{animal_num}{session_num}{region_num,other_region_num}{i}(:,shuffle_num_cell) =  exp(right_all_task1_var{animal_num}{session_num}(sig_coupling_predictor_idx,:)' * all_coeff{animal_num}{session_num}(cell_num,sig_coupling_predictor_idx)')... % Y hat for useful cell predictors
                                        .* exp(B0(cell_num));
                                else
                                    control_cell_right_predictors_predicted_activity{animal_num}{session_num}{region_num,other_region_num}{i}(:,shuffle_num_cell) = nan(size(right_all_task1_var{animal_num}{session_num},2),1);
                                    original_cell_right_predictors_predicted_activity{animal_num}{session_num}{region_num,other_region_num}{i}(:,shuffle_num_cell) = nan(size(right_all_task1_var{animal_num}{session_num},2),1);
                                end
                                
                                % Reconstruct left trials
                                if ~isempty(region_cells) & ~isempty(left_all_task1_var{animal_num}{session_num})
                                    control_cell_left_predictors_predicted_activity{animal_num}{session_num}{region_num,other_region_num}{i}(:,shuffle_num_cell) =  exp(left_all_task1_var{animal_num}{session_num}(sig_outside_coupling_predictor_idx,:)' * all_coeff{animal_num}{session_num}(cell_num,sig_outside_coupling_predictor_idx)')... % Y hat for useful cell predictors
                                        .* exp(B0(cell_num));
                                    original_cell_left_predictors_predicted_activity{animal_num}{session_num}{region_num,other_region_num}{i}(:,shuffle_num_cell) =  exp(left_all_task1_var{animal_num}{session_num}(sig_coupling_predictor_idx,:)' * all_coeff{animal_num}{session_num}(cell_num,sig_coupling_predictor_idx)')... % Y hat for useful cell predictors
                                        .* exp(B0(cell_num));
                                else
                                    control_cell_left_predictors_predicted_activity{animal_num}{session_num}{region_num,other_region_num}{i}(:,shuffle_num_cell) = nan(size(left_all_task1_var{animal_num}{session_num},2),1);
                                    original_cell_left_predictors_predicted_activity{animal_num}{session_num}{region_num,other_region_num}{i}(:,shuffle_num_cell) = nan(size(left_all_task1_var{animal_num}{session_num},2),1);
                                end
                                
                                clear rand_selected_region_task_sig_cell rand_selected_all_region_task_sig_cell rand_selected_non_region_task_sig_cell
                            end
                        end
                        
                        clear control_predictor_idx removed_predictor_idx
                        epoch_duration = length(epoch_start:epoch_end);
                        % Start arranging cells back into its trials
                        num_right_trials = size(control_cell_right_predictors_predicted_activity{animal_num}{session_num}{region_num,other_region_num}{i},1)/epoch_duration;
                        num_left_trials = size(control_cell_left_predictors_predicted_activity{animal_num}{session_num}{region_num,other_region_num}{i},1)/epoch_duration;
                        
                        % Reshape time series into trials
                        trial_control_cell_right_predictors_activity{animal_num}{session_num}{region_num,other_region_num}{i} = reshape(control_cell_right_predictors_predicted_activity{animal_num}{session_num}{region_num,other_region_num}{i},[epoch_duration,num_right_trials,fixed_num_shuffle]);
                        trial_original_cell_right_predictors_activity{animal_num}{session_num}{region_num,other_region_num}{i} = reshape(original_cell_right_predictors_predicted_activity{animal_num}{session_num}{region_num,other_region_num}{i},[epoch_duration,num_right_trials,fixed_num_shuffle]);
                        trial_control_cell_left_predictors_activity{animal_num}{session_num}{region_num,other_region_num}{i} = reshape(control_cell_left_predictors_predicted_activity{animal_num}{session_num}{region_num,other_region_num}{i},[epoch_duration,num_left_trials,fixed_num_shuffle]);
                        trial_original_cell_left_predictors_activity{animal_num}{session_num}{region_num,other_region_num}{i} = reshape(original_cell_left_predictors_predicted_activity{animal_num}{session_num}{region_num,other_region_num}{i},[epoch_duration,num_left_trials,fixed_num_shuffle]);
                        
                        clear region_predictor_idx randomly_shuffled_to_remove_region_predictor_idx removed_predictor_idx control_to_remove_predictor_idx control_predictor_idx...
                            control_cell_left_predictors_predicted_activity original_cell_left_predictors_predicted_activity...
                            control_cell_right_predictors_predicted_activity  original_cell_right_predictors_predicted_activity
                    end
                end
            end
        end
        
    end
    
    % Save reconstructed ablated activity
    ablated_data.reconstructed_activity.trial_control_cell_right_predictors_activity = trial_control_cell_right_predictors_activity;
    ablated_data.reconstructed_activity.trial_original_cell_right_predictors_activity = trial_original_cell_right_predictors_activity;
    ablated_data.reconstructed_activity.trial_control_cell_left_predictors_activity = trial_control_cell_left_predictors_activity;
    ablated_data.reconstructed_activity.trial_original_cell_left_predictors_activity = trial_original_cell_left_predictors_activity;
    
    clear task1_var_temp right_all_task1_var right_all_resp_var_trial_concat left_all_resp_var_trial_concat left_all_task1_var both_all_resp_var_trial_concat both_all_task1_var...
        control_cell_left_predictors_predicted_activity  original_cell_left_predictors_predicted_activity...
        control_cell_right_predictors_predicted_activity   original_cell_right_predictors_predicted_activity
end

%% Part 2: Compute choice mode from reconstructed activity
for animal_num  = 1:14
    session_num = 1;
    for region_num = region_idx
        for other_region_num =region_idx
            
            if ~isempty(ablated_data.reconstructed_activity.trial_control_cell_right_predictors_activity{animal_num}{session_num}{region_num,other_region_num})
                for i = 1:length(find(abs(cellfun(@isempty,ablated_data.reconstructed_activity.trial_control_cell_right_predictors_activity{animal_num}{session_num}{region_num,other_region_num})-1)))
                    cell_idx = find(abs(cellfun(@isempty,ablated_data.reconstructed_activity.trial_control_cell_right_predictors_activity{animal_num}{session_num}{region_num,other_region_num})-1));
                    cell_num = cell_idx(i);
                    % Rearrange cell to get frame x trial x cell
                    reshaped_control_left_predicted_activity{animal_num}{session_num}{region_num,other_region_num}(:,:,:,i) = ablated_data.reconstructed_activity.trial_control_cell_left_predictors_activity{animal_num}{session_num}{region_num,other_region_num}{cell_num};
                    reshaped_original_left_predicted_activity{animal_num}{session_num}{region_num,other_region_num}(:,:,:,i) = ablated_data.reconstructed_activity.trial_original_cell_left_predictors_activity{animal_num}{session_num}{region_num,other_region_num}{cell_num};
                    
                    reshaped_control_right_predicted_activity{animal_num}{session_num}{region_num,other_region_num}(:,:,:,i) = ablated_data.reconstructed_activity.trial_control_cell_right_predictors_activity{animal_num}{session_num}{region_num,other_region_num}{cell_num};
                    reshaped_original_right_predicted_activity{animal_num}{session_num}{region_num,other_region_num}(:,:,:,i) = ablated_data.reconstructed_activity.trial_original_cell_right_predictors_activity{animal_num}{session_num}{region_num,other_region_num}{cell_num};
                end
            else
                reshaped_control_left_predicted_activity{animal_num}{session_num}{region_num,other_region_num}= [];
                reshaped_original_left_predicted_activity{animal_num}{session_num}{region_num,other_region_num} = [];
                reshaped_control_right_predicted_activity{animal_num}{session_num}{region_num,other_region_num} = [];
                reshaped_original_right_predicted_activity{animal_num}{session_num}{region_num,other_region_num} = [];
            end
        end
    end
    
    %% Get choice mode population activity by projecting along choice axis
    region= 8;
    choice_epoch_idx = round(fs_image*3.5):round(fs_image*4);
    session_idx = 1;
    for session_num = session_idx
        for region_num = region_idx
            for other_region_num = region_idx
                % Valid region with cells
                if ~isempty(reshaped_original_left_predicted_activity{animal_num}{session_num}{region_num,other_region_num})
                    if length(find(reshaped_original_left_predicted_activity{animal_num}{session_num}{region_num,other_region_num}(1,1,1,:)>0)) >= 3 % At least 3 cells in a given region.
                        
                        % Get num cells x 1 vector
                        left_original_choice = squeeze(mean(mean(reshaped_original_left_predicted_activity{animal_num}{session_num}{region_num,other_region_num}(choice_epoch_idx,:,:,:),1),2));
                        left_control_choice = squeeze(mean(mean(reshaped_control_left_predicted_activity{animal_num}{session_num}{region_num,other_region_num}(choice_epoch_idx,:,:,:),1),2));
                        
                        left_original_choice = left_original_choice(:,find(left_original_choice(1,:)~=0));
                        left_control_choice = left_control_choice(:,find(left_control_choice(1,:)~=0));
                        
                        left_original_choice = left_original_choice(:,~isnan(left_original_choice(1,:)));
                        left_control_choice = left_control_choice(:,~isnan(left_control_choice(1,:)));
                        
                        % Get num cells x 1 vector
                        right_original_choice = squeeze(mean(mean(reshaped_original_right_predicted_activity{animal_num}{session_num}{region_num,other_region_num}(choice_epoch_idx,:,:,:),1),2));
                        right_control_choice = squeeze(mean(mean(reshaped_control_right_predicted_activity{animal_num}{session_num}{region_num,other_region_num}(choice_epoch_idx,:,:,:),1),2));
                        
                        right_original_choice = right_original_choice(:,find(right_original_choice(1,:)~=0));
                        right_control_choice = right_control_choice(:,find(right_control_choice(1,:)~=0));
                        
                        right_original_choice = right_original_choice(:,~isnan(right_original_choice(1,:)));
                        right_control_choice = right_control_choice(:,~isnan(right_control_choice(1,:)));
                        
                        % Choice diff - cell x 1 (Averaged across frames)
                        original_choice_diff = right_original_choice - left_original_choice;
                        control_choice_diff = right_control_choice - left_control_choice;
                        
                        original_norm_choice_diff = original_choice_diff./((sum(original_choice_diff.^2,2)).^0.5);
                        control_norm_choice_diff = control_choice_diff./((sum(control_choice_diff.^2,2)).^0.5);
                        
                        % Project to choice mode.
                        % All resp - cell x trial x frame x shuffle
                        original_left_all_resp = permute(reshaped_original_left_predicted_activity{animal_num}{session_num}{region_num,other_region_num},[4 2 1 3]);
                        control_left_all_resp = permute(reshaped_control_left_predicted_activity{animal_num}{session_num}{region_num,other_region_num},[4 2 1 3]);
                        
                        original_right_all_resp = permute(reshaped_original_right_predicted_activity{animal_num}{session_num}{region_num,other_region_num},[4 2 1 3]);
                        control_right_all_resp = permute(reshaped_control_right_predicted_activity{animal_num}{session_num}{region_num,other_region_num},[4 2 1 3]);
                        
                        original_left_all_resp = original_left_all_resp(find(original_left_all_resp(:,1,1,1) ~=0),:,:,:);
                        control_left_all_resp = control_left_all_resp(find(control_left_all_resp(:,1,1,1) ~=0),:,:,:);
                        
                        original_right_all_resp = original_right_all_resp(find(original_right_all_resp(:,1,1,1) ~=0),:,:,:);
                        control_right_all_resp = control_right_all_resp(find(control_right_all_resp(:,1,1,1) ~=0),:,:,:);
                        
                        original_left_all_resp = original_left_all_resp(~isnan(original_left_all_resp(:,1,1,1)),:,:,:);
                        control_left_all_resp = control_left_all_resp(~isnan(control_left_all_resp(:,1,1,1)),:,:,:);
                        
                        original_right_all_resp = original_right_all_resp(~isnan(original_right_all_resp(:,1,1,1)),:,:,:);
                        control_right_all_resp = control_right_all_resp(~isnan(control_right_all_resp(:,1,1,1)),:,:,:);
                        
                        for trial_num = 1:size(original_left_all_resp,2)
                            % trials x frame - for each animal
                            for shuffle_num = 1:size(original_left_all_resp,4)
                                original_left_choice_mode{animal_num}{session_num}{region_num,other_region_num}(trial_num,:,shuffle_num) = (squeeze(original_left_all_resp(:,trial_num,:,shuffle_num)))'*original_norm_choice_diff(shuffle_num,:)';
                                control_left_choice_mode{animal_num}{session_num}{region_num,other_region_num}(trial_num,:,shuffle_num) = (squeeze(control_left_all_resp(:,trial_num,:,shuffle_num)))'*control_norm_choice_diff(shuffle_num,:)';
                            end
                        end
                        for trial_num = 1:size(original_right_all_resp,2)
                            for shuffle_num = 1:size(original_right_all_resp,4)
                                original_right_choice_mode{animal_num}{session_num}{region_num,other_region_num}(trial_num,:,shuffle_num) = (squeeze(original_right_all_resp(:,trial_num,:,shuffle_num)))'*original_norm_choice_diff(shuffle_num,:)';
                                control_right_choice_mode{animal_num}{session_num}{region_num,other_region_num}(trial_num,:,shuffle_num) = (squeeze(control_right_all_resp(:,trial_num,:,shuffle_num)))'*control_norm_choice_diff(shuffle_num,:)';
                            end
                        end
                    else
                        original_left_choice_mode{animal_num}{session_num}{region_num,other_region_num} = [];%nan(size(original_left_all_resp,2),size(original_left_all_resp,3));
                        control_left_choice_mode{animal_num}{session_num}{region_num,other_region_num} = [];%nan(size(control_left_all_resp,2),size(original_left_all_resp,3));
                        original_right_choice_mode{animal_num}{session_num}{region_num,other_region_num} = [];%nan(size(original_right_all_resp,2),size(original_left_all_resp,3));
                    end
                    clear choice_diff norm_choice_diff left_original_choice left_control_choice right_control_choice all_resp norm_choice_diff
                else
                    original_left_choice_mode{animal_num}{session_num}{region_num,other_region_num} = [];%nan(size(original_left_all_resp,2),size(original_left_all_resp,3));
                    control_left_choice_mode{animal_num}{session_num}{region_num,other_region_num} = [];%nan(size(control_left_all_resp,2),size(original_left_all_resp,3));
                    
                    original_right_choice_mode{animal_num}{session_num}{region_num,other_region_num} = [];%nan(size(original_right_all_resp,2),size(original_left_all_resp,3));
                    control_right_choice_mode{animal_num}{session_num}{region_num,other_region_num} = [];%nan(size(control_left_all_resp,2),size(original_left_all_resp,3));
                end
            end
        end
    end
    
    % Concatenate all the choice mode population activity across all animals
    for region_num = region_idx
        for other_region_num = region_idx
            concat_original_right_choice_mode{session_num}{region_num,other_region_num}{animal_num}  = [];
            concat_control_right_choice_mode{session_num}{region_num,other_region_num}{animal_num}  = [];
            
            concat_original_left_choice_mode{session_num}{region_num,other_region_num}{animal_num}  = [];
            concat_control_left_choice_mode{session_num}{region_num,other_region_num}{animal_num}  = [];
            
            if ~isempty(original_right_choice_mode{animal_num}{session_num}{region_num,other_region_num})
                concat_original_right_choice_mode{session_num}{region_num,other_region_num}{animal_num} = original_right_choice_mode{animal_num}{session_num}{region_num,other_region_num};
                concat_control_right_choice_mode{session_num}{region_num,other_region_num}{animal_num} = control_right_choice_mode{animal_num}{session_num}{region_num,other_region_num};
                
                concat_original_left_choice_mode{session_num}{region_num,other_region_num}{animal_num} = original_left_choice_mode{animal_num}{session_num}{region_num,other_region_num};
                concat_control_left_choice_mode{session_num}{region_num,other_region_num}{animal_num} = control_left_choice_mode{animal_num}{session_num}{region_num,other_region_num};
            end
        end
    end
end
% save data in struct
ablated_data.choice_mode.concat_original_right_choice_mode = concat_original_right_choice_mode;
ablated_data.choice_mode.concat_original_right_choice_mode = concat_control_right_choice_mode;
ablated_data.choice_mode.concat_original_right_choice_mode = concat_original_left_choice_mode;
ablated_data.choice_mode.concat_original_right_choice_mode = concat_control_left_choice_mode;
end