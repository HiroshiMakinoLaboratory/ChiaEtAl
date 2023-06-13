%% Analyse retained vs lost coupling by selectively reconstructing activity
function retained_elim_data = get_retained_eliminated_activity(parameters,data)
% Initialise variables
epoch_start = round(parameters.fs_image*3); % +3s from start of (-4s) ITI
epoch_end = round(parameters.fs_image*8); % 1s ITI + 1s stimulus + 2s delay + 1s action
animal_idx = 1:13;
fixed_num_shuffle = 100; % Subsampling of functional coupling is performed 100 times

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
all_resp_var_trial_concat = data.activity.all_resp_var_trial_concat; % Neural activity
all_trial_choice = data.trial_idx.all_trial_choice;
all_valid_cells = data.cell_idx.all_valid_cells;
all_region_labels = data.cell_idx.all_region_labels;
all_predictor_labels = data.cell_idx.all_predictor_labels;
coupling_sig_cells = data.cell_idx.coupling_sig_cells;
all_coeff = data.cell_idx.all_coeff;
all_GLM = data.all_GLM;
target_region_num = [1:8]; % Use all regions

for session_num = [2] % Reconstruct intermediate session
    % Get reconstructed predicted activity based on coupling predictors
    for animal_num = animal_idx
        
        clear task1_var_temp left_trial_idx right_trial_idx
        % Get Task1_var_temp which are the normalised predictors (i.e.task + coupling predictors)
        % Using 1s ITI + 1s stimulus + 2s delay + 1s action
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
        
        %% Splitting trial types (correct left, incorrect left, etc)
        % Get trial types based on animal choice regardless of stimulus
        left_trial_idx = sort([find(all_trial_choice{animal_num}{session_num} == 1)])';
        right_trial_idx = sort([find(all_trial_choice{animal_num}{session_num} == 2)])';
        
        clear left_all_resp_var_trial_concat left_all_task1_var right_all_resp_var_trial_concat right_all_task1_var
        
        % Both choices
        both_all_task1_var{animal_num}{session_num} = []; % Concat predicted activity
        for trial_num = 1:size(all_resp_var_trial_concat{animal_num}{session_num},2)
            both_all_task1_var{animal_num}{session_num}  = [both_all_task1_var{animal_num}{session_num} task1_var_temp{animal_num}{session_num}{trial_num}'];
        end
        
        %% Get predicted activity using coeff
        % Inititalise variables
        B0 = all_GLM{animal_num}{session_num}.B0;
        for i = 1: length(all_region_labels{animal_num}{session_num})
            retained_naive_predictors_predicted_activity{animal_num}{session_num}{i} = [];
            lost_naive_predictors_predicted_activity{animal_num}{session_num}{i} = [];
        end
        
        % Get predicted activity by reconstructing the retained or eliminated activity
        % Reconstruct control vs silenced population
        if length(find(ismember(all_region_labels{animal_num}{session_num},target_region_num))) ~= 0
            for i = 1:length( find(ismember(all_region_labels{animal_num}{session_num},target_region_num)) )
                cell_num = i;
                % Get animal specific predictor labels
                predictor_labels = all_predictor_labels{animal_num}{session_num};
                predictor_idx = find(~isnan(predictor_labels));
                
                if ~isempty(all_valid_cells{animal_num}{session_num})
                    clear size_combi_vector
                    rng(2020);
                    for shuffle_num_cell = 1:fixed_num_shuffle % This is to randomly select cells to be removed
                        
                        % Intermediate session
                        naive_sig_coupling = coupling_sig_cells{animal_num}{session_num}{cell_num};
                        
                        % Expert session
                        expert_sig_coupling = coupling_sig_cells{animal_num}{3}{cell_num};
                        
                        % Coupling predictors found in both intermediate and expert
                        retained_naive_coupling = intersect(intersect(naive_sig_coupling,expert_sig_coupling),all_valid_cells{animal_num}{session_num});
                        
                        % Coupling predictors found in intermediate but NOT expert
                        lost_naive_coupling = intersect(naive_sig_coupling(~ismember(naive_sig_coupling,retained_naive_coupling)),all_valid_cells{animal_num}{session_num});
                        
                        total_num_retained_naive_coupling{animal_num}{session_num}(i) = length(retained_naive_coupling);
                        total_num_lost_naive_coupling{animal_num}{session_num}(i) = length(lost_naive_coupling);
                        
                        % Randomly shuffle to ensure the same number of coupling
                        min_num_coupling = min([length(retained_naive_coupling) length(lost_naive_coupling)]);
                        
                        % Get predictor index (because each cell coupling has 2 predictors t-1 and t-2)
                        retained_naive_coupling_predictor_idx = find(ismember(predictor_labels ,retained_naive_coupling(randsample(length(retained_naive_coupling),min_num_coupling)) ));
                        lost_naive_coupling_coupling_predictor_idx = find(ismember(predictor_labels ,lost_naive_coupling(randsample(length(lost_naive_coupling),min_num_coupling)) ));
                        
                        % Reconstruct using selected cell-coupling predictors
                        if ~isempty(both_all_task1_var{animal_num}{session_num})
                            retained_naive_predictors_predicted_activity{animal_num}{session_num}{i}(:,shuffle_num_cell) =  exp(both_all_task1_var{animal_num}{session_num}(retained_naive_coupling_predictor_idx,:)' * all_coeff{animal_num}{session_num}(cell_num,retained_naive_coupling_predictor_idx)')... % Y hat for useful cell predictors
                                .* exp(B0(cell_num));
                            lost_naive_predictors_predicted_activity{animal_num}{session_num}{i}(:,shuffle_num_cell) =  exp(both_all_task1_var{animal_num}{session_num}(lost_naive_coupling_coupling_predictor_idx,:)' * all_coeff{animal_num}{session_num}(cell_num,lost_naive_coupling_coupling_predictor_idx)')... % Y hat for useful cell predictors
                                .* exp(B0(cell_num));
                        else
                            retained_naive_predictors_predicted_activity{animal_num}{session_num}{i}(:,shuffle_num_cell) = nan(size(both_all_task1_var{animal_num}{session_num},2),1);
                            lost_naive_predictors_predicted_activity{animal_num}{session_num}{i}(:,shuffle_num_cell) = nan(size(both_all_task1_var{animal_num}{session_num},2),1);
                        end
                        clear rand_selected_region_task_sig_cell   rand_selected_all_region_task_sig_cell  rand_selected_non_region_task_sig_cell
                    end
                    
                    % reshape time series back into trials
                    for hide = 1
                        clear control_predictor_idx removed_predictor_idx
                        epoch_duration = length(epoch_start:epoch_end);
                        
                        % Start arranging cells back into its trials
                        num_both_trials = size(retained_naive_predictors_predicted_activity{animal_num}{session_num}{i},1)/epoch_duration;
                        
                        % Reshape time series into trials
                        % Frames x shuffle
                        trial_retained_naive_predictors_predicted_activity{animal_num}{i} = reshape(retained_naive_predictors_predicted_activity{animal_num}{session_num}{i},[epoch_duration,num_both_trials,fixed_num_shuffle]);
                        trial_lost_naive_predictors_predicted_activity{animal_num}{i} = reshape(lost_naive_predictors_predicted_activity{animal_num}{session_num}{i},[epoch_duration,num_both_trials,fixed_num_shuffle]);
                        
                        % Check reshape is correct
                        a = trial_retained_naive_predictors_predicted_activity{animal_num}{i}(:,2,2);
                        b = retained_naive_predictors_predicted_activity{animal_num}{session_num}{i}( (epoch_duration+1 : epoch_duration*2),2 );
                        assert(isequal(a,b));
                        
                        clear region_predictor_idx randomly_shuffled_to_remove_region_predictor_idx removed_predictor_idx control_to_remove_predictor_idx control_predictor_idx
                    end
                    clear control_cell_both_predictors_predicted_activity  silenced_cell_both_predictors_predicted_activity  original_cell_both_predictors_predicted_activity...
                        control_cell_left_predictors_predicted_activity silenced_cell_left_predictors_predicted_activity original_cell_left_predictors_predicted_activity...
                        control_cell_right_predictors_predicted_activity  silenced_cell_right_predictors_predicted_activity original_cell_right_predictors_predicted_activity...
                        lost_interm_predictors_predicted_activity  retained_naive_predictors_predicted_activity retained_interm_predictors_predicted_activity lost_naive_predictors_predicted_activity
                end
            end
            
            % Reshape into:
            % cell x trial x frame x shuffle
            for i = 1:size(trial_retained_naive_predictors_predicted_activity{animal_num},2)
                reshaped_retained_naive_predicted_activity{animal_num}{session_num}(i,:,:,:) = permute(trial_retained_naive_predictors_predicted_activity{animal_num}{i},[2,1,3]);
                reshaped_lost_naive_predicted_activity{animal_num}{session_num}(i,:,:,:) = permute(trial_lost_naive_predictors_predicted_activity{animal_num}{i},[2,1,3]);
            end
            
        else
            reshaped_retained_naive_predicted_activity{animal_num}{session_num} = [];
            reshaped_lost_naive_predicted_activity{animal_num}{session_num} = [];
        end
    end
end

% Save variables intro struct
retained_elim_data.reconstructed_activity.reshaped_retained_naive_predicted_activity = reshaped_retained_naive_predicted_activity;
retained_elim_data.reconstructed_activity.reshaped_lost_naive_predicted_activity = reshaped_lost_naive_predicted_activity;
retained_elim_data.num_cells.total_num_lost_naive_coupling = total_num_lost_naive_coupling;
retained_elim_data.num_cells.total_num_retained_naive_coupling = total_num_retained_naive_coupling;
end