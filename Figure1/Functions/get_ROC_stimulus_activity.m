%% Get ROC activity without averaging across session, time and cells.
function ROC_data = get_ROC_stimulus_activity(data,stimulus_or_choice)

if strcmp(stimulus_or_choice,'stimulus')
    raw_left_trial_activity = data.activity.raw_left_trial_activity;
    raw_right_trial_activity = data.activity.raw_right_trial_activity;
elseif strcmp(stimulus_or_choice,'choice')
    raw_left_trial_activity = data.activity.raw_left_choice_activity;
    raw_right_trial_activity = data.activity.raw_right_choice_activity;
end

clear right_stim_resp left_stim_resp
for animal_num = 1:size(raw_left_trial_activity,2)
    for session_num = 1:size(raw_left_trial_activity{animal_num},2)
        for region_num = 1:8
            region_idx = intersect(find(data.cell_idx.all_region_labels{animal_num}{session_num} == region_num),find(data.cell_idx.all_valid_cells{animal_num}{session_num}));
            
            if ~isempty(raw_right_trial_activity{animal_num}{session_num})
                right_stim_resp{animal_num}{session_num}{region_num} = permute(raw_right_trial_activity{animal_num}{session_num}(:,:,region_idx),[3 1 2]);
            else
                right_stim_resp{animal_num}{session_num}{region_num} =  [];
            end
            
            if ~isempty(raw_left_trial_activity{animal_num}{session_num})
                left_stim_resp{animal_num}{session_num}{region_num} = permute(raw_left_trial_activity{animal_num}{session_num}(:,:,region_idx),[3 1 2]);
            else
                left_stim_resp{animal_num}{session_num}{region_num}=  [];
            end
        end
    end
end

num_shuffles = 100;
clear cell_auc  shuffled_cell_auc
for animal_num = 1:13
    for session_num = 1:3
        for region_num = 1:8
            % cells x trials x time series
            if size(left_stim_resp{animal_num}{session_num}{region_num},1) > 5
                for cell_num = 1:size(left_stim_resp{animal_num}{session_num}{region_num},1)
                    if (size(right_stim_resp{animal_num}{session_num}{region_num},2) >= 5 & (size(left_stim_resp{animal_num}{session_num}{region_num},2) >= 5))
                        for time_point = 1:size(right_stim_resp{animal_num}{session_num}{region_num}(cell_num,:,:),3)
                            % Each datapoint is a trial
                            trial_a = right_stim_resp{animal_num}{session_num}{region_num}(cell_num,:,time_point);
                            trial_b = left_stim_resp{animal_num}{session_num}{region_num}(cell_num,:,time_point);
                            
                            roc_data = roc_curve(trial_a,trial_b);
                            cell_auc{animal_num}{session_num}{region_num}(cell_num,time_point) = abs(roc_data.param.AROC-0.5) .*2;
                            
                            % Get null distribution
                            rng(2020)
                            for shuffle_num = 1:num_shuffles
                                % shuffle the labels
                                concat_trials = [trial_a trial_b];
                                
                                rand_trial_a_idx = randsample(1:(length(trial_a)+length(trial_b)),length(trial_a));
                                rand_trial_b_idx = randsample(1:(length(trial_a)+length(trial_b)),length(trial_b));
                                
                                rand_trial_a = concat_trials(rand_trial_a_idx);
                                rand_trial_b = concat_trials(rand_trial_b_idx);
                                
                                roc_data = roc_curve(rand_trial_a,rand_trial_b);
                                shuffled_cell_auc{animal_num}{session_num}{region_num}(cell_num,time_point,shuffle_num) = abs(roc_data.param.AROC-0.5) .*2;
                            end
                        end
                    else
                        cell_auc{animal_num}{session_num}{region_num} = [];
                        shuffled_cell_auc{animal_num}{session_num}{region_num} = [];
                    end
                end
            else
                cell_auc{animal_num}{session_num}{region_num} =[];
                shuffled_cell_auc{animal_num}{session_num}{region_num} = [];
            end
        end
    end
end

ROC_data.cell_auc = cell_auc;
ROC_data.shuffled_cell_auc = shuffled_cell_auc;

%% Concat across animals
data_analysis = cell_auc;
shuffled_data_analysis = shuffled_cell_auc;

for session_num = 1:3
    for region_num = 1:8
        concat_ROC_data{session_num}{region_num} = [];
        concat_shuffled_ROC_data{session_num}{region_num} = [];
        
        for animal_num = 1:13
            concat_ROC_data{session_num}{region_num} = [concat_ROC_data{session_num}{region_num} ; data_analysis{animal_num}{session_num}{region_num}];
            concat_shuffled_ROC_data{session_num}{region_num} = [concat_shuffled_ROC_data{session_num}{region_num} ; shuffled_data_analysis{animal_num}{session_num}{region_num}];
        end
    end
end

ROC_data.concat_auc = concat_ROC_data;
ROC_data.concat_shuffled_auc_data = concat_shuffled_ROC_data;
end