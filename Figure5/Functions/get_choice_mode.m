%% Get choice mode from retained or eliminated reconstructed activity
function retained_elim_data = get_choice_mode(parameters,data,retained_elim_data)
% Load variables
reshaped_retained_naive_predicted_activity = retained_elim_data.reconstructed_activity.reshaped_retained_naive_predicted_activity;
reshaped_lost_naive_predicted_activity = retained_elim_data.reconstructed_activity.reshaped_lost_naive_predicted_activity;

all_trial_choice = data.trial_idx.all_trial_choice;
delay_right_sig_cells = data.cell_idx.delay_right_sig_cells;
fs_image = parameters.fs_image;
choice_epoch_idx = round(fs_image*3.5):round(fs_image*4);
trial_threshold = 4; % Each region requires this number of cells in order for it to be valid

for session_num = 2 % Intermediate session
    for animal_num = 1:13
        if ~isempty(reshaped_retained_naive_predicted_activity{animal_num}{session_num})
            
            % Cell x trial x frame x shuffle
            retained_naive_left_choice_mode{animal_num}{session_num}= [];
            lost_naive_left_choice_mode{animal_num}{session_num}= [];
            
            retained_naive_right_choice_mode{animal_num}{session_num}= [];
            lost_naive_right_choice_mode{animal_num}{session_num}= [];
            
            retained_naive_both_choice_mode{animal_num}{session_num}= [];
            lost_naive_both_choice_mode{animal_num}{session_num}= [];
            
            left_choice_idx = find(all_trial_choice{animal_num}{session_num} == 1);
            right_choice_idx = find(all_trial_choice{animal_num}{session_num} == 2);
            
            % Use delay right encoding cells
            task_variable_cell = delay_right_sig_cells{animal_num}{session_num};
            
            if length(left_choice_idx) > trial_threshold & length(right_choice_idx) > trial_threshold & length(task_variable_cell) >2
                
                % Get num cells x 1 vector
                % Average cell x trial x frame x shuffle to remove frame and tiral
                left_retained_naive_choice = squeeze(mean(mean(reshaped_retained_naive_predicted_activity{animal_num}{session_num}(task_variable_cell,left_choice_idx,choice_epoch_idx,:),2),3));
                left_lost_naive_choice = squeeze(mean(mean(reshaped_lost_naive_predicted_activity{animal_num}{session_num}(task_variable_cell,left_choice_idx,choice_epoch_idx,:),2),3));
                
                % Get num cells x 1 vector
                right_retained_naive_choice = squeeze(mean(mean(reshaped_retained_naive_predicted_activity{animal_num}{session_num}(task_variable_cell,right_choice_idx,choice_epoch_idx,:),2),3));
                right_lost_naive_choice = squeeze(mean(mean(reshaped_lost_naive_predicted_activity{animal_num}{session_num}(task_variable_cell,right_choice_idx,choice_epoch_idx,:),2),3));
                
                % Choice diff - cell x 1 (Averaged across frames)
                retained_naive_choice_diff = right_retained_naive_choice - left_retained_naive_choice;
                lost_naive_choice_diff = right_lost_naive_choice - left_lost_naive_choice;
                
                % L2 norm
                retained_naive_norm_choice_diff = retained_naive_choice_diff./((sum(retained_naive_choice_diff.^2,1)).^0.5);
                lost_naive_norm_choice_diff = lost_naive_choice_diff./((sum(lost_naive_choice_diff.^2,1)).^0.5);
                
                % Project to choice mode.
                % All resp - cell x trial x frame x shuffle
                retained_naive_all_resp = reshaped_retained_naive_predicted_activity{animal_num}{session_num}(task_variable_cell,:,:,:);
                lost_naive_all_resp = reshaped_lost_naive_predicted_activity{animal_num}{session_num}(task_variable_cell,:,:,:);
                
                for i = 1:length(left_choice_idx)
                    trial_num = left_choice_idx(i);
                    % trials x frame - for each animal
                    for shuffle_num = 1:size(reshaped_retained_naive_predicted_activity{animal_num}{session_num},4)
                        retained_left_choice_mode{animal_num}{session_num}(i,:,shuffle_num) = (squeeze(retained_naive_all_resp(:,trial_num,:,shuffle_num)))'*retained_naive_norm_choice_diff(:,shuffle_num);
                        lost_left_choice_mode{animal_num}{session_num}(i,:,shuffle_num) = (squeeze(lost_naive_all_resp(:,trial_num,:,shuffle_num)))'*lost_naive_norm_choice_diff(:,shuffle_num);
                    end
                end
                
                for i = 1:length(right_choice_idx)
                    trial_num = right_choice_idx(i);
                    % trials x frame - for each animal
                    for shuffle_num = 1:size(reshaped_retained_naive_predicted_activity{animal_num}{session_num},4)
                        retained_right_choice_mode{animal_num}{session_num}(i,:,shuffle_num) = (squeeze(retained_naive_all_resp(:,trial_num,:,shuffle_num)))'*retained_naive_norm_choice_diff(:,shuffle_num);
                        lost_right_choice_mode{animal_num}{session_num}(i,:,shuffle_num) = (squeeze(lost_naive_all_resp(:,trial_num,:,shuffle_num)))'*lost_naive_norm_choice_diff(:,shuffle_num);
                    end
                end
                
                both_choice_idx = sort(union(right_choice_idx,left_choice_idx));
                for i = 1:180
                    % trials x frame - for each animal
                    for shuffle_num = 1:size(reshaped_retained_naive_predicted_activity{animal_num}{session_num},4)
                        retained_both_choice_mode{animal_num}{session_num}(i,:,shuffle_num) = (squeeze(retained_naive_all_resp(:,trial_num,:,shuffle_num)))'*retained_naive_norm_choice_diff(:,shuffle_num);
                        lost_both_choice_mode{animal_num}{session_num}(i,:,shuffle_num) = (squeeze(lost_naive_all_resp(:,trial_num,:,shuffle_num)))'*lost_naive_norm_choice_diff(:,shuffle_num);
                    end
                end
                
            else
                retained_right_choice_mode{animal_num}{session_num} = [];
                lost_right_choice_mode{animal_num}{session_num} = [];
                retained_left_choice_mode{animal_num}{session_num} =[];
                lost_left_choice_mode{animal_num}{session_num} = [];
                retained_both_choice_mode{animal_num}{session_num}= [];
                lost_both_choice_mode{animal_num}{session_num} = [];
            end
            
        else
            retained_right_choice_mode{animal_num}{session_num} = [];
            lost_right_choice_mode{animal_num}{session_num} = [];
            retained_left_choice_mode{animal_num}{session_num} =[];
            lost_left_choice_mode{animal_num}{session_num} = [];
            retained_both_choice_mode{animal_num}{session_num}= [];
            lost_both_choice_mode{animal_num}{session_num} = [];
        end
        
        clear choice_diff norm_choice_diff left_original_choice left_silenced_choice left_control_choice both_original_choice right_silenced_choice right_control_choice all_resp norm_choice_diff
    end
end

retained_elim_data.choice_mode.retained_right_choice_mode = retained_right_choice_mode;
retained_elim_data.choice_mode.lost_right_choice_mode = lost_right_choice_mode;
retained_elim_data.choice_mode.retained_left_choice_mode = retained_left_choice_mode;
retained_elim_data.choice_mode.lost_left_choice_mode = lost_left_choice_mode;
retained_elim_data.choice_mode.retained_both_choice_mode = retained_both_choice_mode;
retained_elim_data.choice_mode.lost_both_choice_mode = lost_both_choice_mode;
end