%% Get ROC activity without averaging across session, time and cells.
function ROC_population_data = get_ROC_choice_stim_mode(data,stimulus_or_choice)
if strcmp(stimulus_or_choice,'choice')
    population_activity = data.all_choice_mode_GS;
    right_trial_type = data.choice_right;
    left_trial_type = data.choice_left;
elseif strcmp(stimulus_or_choice,'stimulus')
    population_activity = data.all_sensory_mode_GS;
    right_trial_type = data.trial_num_left;
    left_trial_type = data.trial_num_right;
end
num_shuffle = 100;

clear cell_auc  shuffled_cell_auc
for animal_num = 1:length(population_activity)
    for session_num = 1:length(population_activity{animal_num})
        for region_num = 1:8
            right_trial_idx = right_trial_type{animal_num}{session_num};
            left_trial_idx = left_trial_type{animal_num}{session_num};
            if (~isnan(population_activity{animal_num}{session_num}{region_num}(1,1)) & length(right_trial_idx) >= 5 & length(left_trial_idx) >= 5 )
                % cells x trials x time series
                for time_point = 1:size(population_activity{animal_num}{session_num}{region_num},2)
                    % Each datapoint is a trial
                    trial_a = population_activity{animal_num}{session_num}{region_num}(right_trial_idx,time_point);
                    trial_b = population_activity{animal_num}{session_num}{region_num}(left_trial_idx,time_point);
                    
                    roc_data = roc_curve(trial_a,trial_b);
                    cell_auc{animal_num}{session_num}{region_num}(time_point) = abs(roc_data.param.AROC-0.5) .*2;
                    
                    % Get null distribution
                    rng(2020)
                    for shuffle_num = 1:num_shuffle
                        % shuffle the labels
                        concat_trials = [trial_a' trial_b'];
                        
                        rand_trial_a_idx = randsample(1:(length(trial_a)+length(trial_b)),length(trial_a));
                        rand_trial_b_idx = randsample(1:(length(trial_a)+length(trial_b)),length(trial_b));
                        
                        rand_trial_a = concat_trials(rand_trial_a_idx);
                        rand_trial_b = concat_trials(rand_trial_b_idx);
                        
                        roc_data = roc_curve(rand_trial_a,rand_trial_b);
                        shuffled_cell_auc{animal_num}{session_num}{region_num}(shuffle_num,time_point) = abs(roc_data.param.AROC-0.5) .*2;
                    end
                end
            else
                cell_auc{animal_num}{session_num}{region_num} = nan(1,size(population_activity{1}{1}{2},2) );
                shuffled_cell_auc{animal_num}{session_num}{region_num} = nan(num_shuffle,size(population_activity{1}{1}{2},2) );
            end
        end
    end
end

ROC_population_data.cell_auc = cell_auc;
ROC_population_data.shuffled_cell_auc = shuffled_cell_auc;

%% Concat across animals

data_analysis = ROC_population_data.cell_auc;
shuffled_data_analysis = ROC_population_data.shuffled_cell_auc;
clear concat_ROC_data  concat_shuffled_ROC_data
for region_num = 1:8
    
    concat_ROC_data{1}{region_num} = [];
    concat_shuffled_ROC_data{1}{region_num} = [];
    % Naive
    session_idx = [1 2];
    for animal_num = 1:length(population_activity)
        concat_ROC_data{1}{region_num} = [concat_ROC_data{1}{region_num} ; data_analysis{animal_num}{session_idx(1)}{region_num} ; data_analysis{animal_num}{session_idx(2)}{region_num}];
        concat_shuffled_ROC_data{1}{region_num} = [concat_shuffled_ROC_data{1}{region_num} ; nanmean(shuffled_data_analysis{animal_num}{session_idx(1)}{region_num},1) ; nanmean(shuffled_data_analysis{animal_num}{session_idx(2)}{region_num},1)];
    end
    
    concat_ROC_data{2}{region_num} = [];
    concat_shuffled_ROC_data{2}{region_num} = [];
    % Intermediate
    session_idx = [3 4];
    for animal_num = 1:length(population_activity)
        concat_ROC_data{2}{region_num} = [concat_ROC_data{2}{region_num} ; data_analysis{animal_num}{session_idx(1)}{region_num} ; data_analysis{animal_num}{session_idx(2)}{region_num}];
        concat_shuffled_ROC_data{2}{region_num} = [concat_shuffled_ROC_data{2}{region_num} ; nanmean(shuffled_data_analysis{animal_num}{session_idx(1)}{region_num},1) ; nanmean(shuffled_data_analysis{animal_num}{session_idx(2)}{region_num},1)];
    end
    
    
    concat_ROC_data{3}{region_num} = [];
    concat_shuffled_ROC_data{3}{region_num} = [];
    % Expert
    session_idx = [5 6];
    for animal_num = 1:length(population_activity)
        concat_ROC_data{3}{region_num} = [concat_ROC_data{3}{region_num} ; data_analysis{animal_num}{session_idx(1)}{region_num} ; data_analysis{animal_num}{session_idx(2)}{region_num}];
        concat_shuffled_ROC_data{3}{region_num} = [concat_shuffled_ROC_data{3}{region_num} ; nanmean(shuffled_data_analysis{animal_num}{session_idx(1)}{region_num},1) ; nanmean(shuffled_data_analysis{animal_num}{session_idx(2)}{region_num},1)];
    end
end
ROC_population_data.concat_auc = concat_ROC_data;
ROC_population_data.concat_shuffled_auc_data = concat_shuffled_ROC_data;
end