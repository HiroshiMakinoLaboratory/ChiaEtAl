%% Get decoding accuracy of retained or eliminated activity
function decoding_accuracy = get_decoding_accuracy(parameters,retained_elim_data)
session_idx = [2]; % intermediate session
trial_threshold = 4; % Number of minimum trials required
fs_image = parameters.fs_image;
choice_epoch_idx = parameters.choice_epoch;

all_trial_choice = retained_elim_data.trials.all_trial_choice;
all_region_labels = retained_elim_data.cell_idx.all_region_labels;
total_num_animals = size(retained_elim_data.trials.all_trial_choice,2);

% Load variables
reshaped_retained_naive_predicted_activity = retained_elim_data.reconstructed_activity.reshaped_retained_naive_predicted_activity;
reshaped_lost_naive_predicted_activity = retained_elim_data.reconstructed_activity.reshaped_lost_naive_predicted_activity;
for animal_num = 1:total_num_animals
    for session_num = session_idx
        % Choice regardless of correct or wrong
        left_choice_idx = find(all_trial_choice{animal_num}{session_num} == 1);
        right_choice_idx = find(all_trial_choice{animal_num}{session_num} == 2);
        
        if length(left_choice_idx) > trial_threshold & length(right_choice_idx) > trial_threshold
            
            % Averaged retained activity
            left_retained_trials = squeeze(mean(reshaped_retained_naive_predicted_activity{animal_num}{session_num}(:,left_choice_idx,choice_epoch_idx,:),3));
            right_retained_trials = squeeze(mean(reshaped_retained_naive_predicted_activity{animal_num}{session_num}(:,right_choice_idx,choice_epoch_idx,:),3));
            
            % Average eliminated activity
            left_lost_trials = squeeze(mean(reshaped_lost_naive_predicted_activity{animal_num}{session_num}(:,left_choice_idx,choice_epoch_idx,:),3));
            right_lost_trials = squeeze(mean(reshaped_lost_naive_predicted_activity{animal_num}{session_num}(:,right_choice_idx,choice_epoch_idx,:),3));
            
            % Each shuffle is 1 subsample for each cell
            for shuffle_num = 1:size(left_lost_trials,3)
                for cell_num = 1: size(left_lost_trials,1)
                    a = left_retained_trials(cell_num,:,shuffle_num);
                    b = right_retained_trials(cell_num,:,shuffle_num);
                    c = left_lost_trials(cell_num,:,shuffle_num);
                    d = right_lost_trials(cell_num,:,shuffle_num);
                    
                    if ~isequal(mean(a),mean(b))
                        auc_params = roc_curve(a,b);
                        auc_retained{session_num}{animal_num}(cell_num,shuffle_num) = auc_params.param.AROC;
                        auc_params = roc_curve(c,d);
                        auc_lost{session_num}{animal_num}(cell_num,shuffle_num) = auc_params.param.AROC;
                    else
                        auc_retained{session_num}{animal_num}(cell_num,shuffle_num) = nan;
                        auc_lost{session_num}{animal_num}(cell_num,shuffle_num) = nan;
                    end
                end
            end
        else
            auc_retained{session_num}{animal_num} = [];
            auc_lost{session_num}{animal_num} = [];
        end
    end
end

% Save decoding accuracy in struct (all regions)
decoding_accuracy.auc_retained = auc_retained;
decoding_accuracy.auc_lost = auc_lost;

%% After performing AUC on all cells across all regions
% Select right choice encoding cells and organise into specific regions
delay_right_sig_cells = retained_elim_data.cell_idx.delay_right_sig_cells;
for session_num = 2
    concat_retained_auc{session_num} = [];
    concat_lost_auc{session_num} = [];
    for region_num = 1:8
        region_concat_retained_auc{session_num}{region_num}  = [];
        region_concat_lost_auc{session_num}{region_num}  = [];
        
    end
    for animal_num = 1:size(decoding_accuracy.auc_retained{session_num},2)
        if ~isempty(decoding_accuracy.auc_retained{session_num}{animal_num})
            % Concat across all sessions
            concat_retained_auc{session_num} = [concat_retained_auc{session_num} mean(decoding_accuracy.auc_retained{session_num}{animal_num},2)'];
            concat_lost_auc{session_num} = [concat_lost_auc{session_num} mean(decoding_accuracy.auc_lost{session_num}{animal_num},2)'];
            
            % Split into regions and average across 100 shuffles for each cell
            for region_num = 1:8
                region_idx_cell = intersect(find(all_region_labels{animal_num}{session_num}==region_num),delay_right_sig_cells{animal_num}{session_num});
                region_concat_retained_auc{session_num}{region_num} = [region_concat_retained_auc{session_num}{region_num} nanmean(decoding_accuracy.auc_retained{session_num}{animal_num}(region_idx_cell,:),2)'];
                region_concat_lost_auc{session_num}{region_num} = [region_concat_lost_auc{session_num}{region_num} nanmean(decoding_accuracy.auc_lost{session_num}{animal_num}(region_idx_cell,:),2)'];
            end
        end
    end
end

% Save decoding accuracy in struct (all regions)
decoding_accuracy.region_concat_retained_auc = region_concat_retained_auc;
decoding_accuracy.region_concat_lost_auc = region_concat_lost_auc;
end