% Get Task-encoding neurons
function task_encoding = get_DREADD_task_encoding_neurons(parameters,DREADD_activity)

cno_all_resp = DREADD_activity.cno.activity.all_resp;
saline_all_resp = DREADD_activity.saline.activity.all_resp;
choice_epoch = parameters.choice_epoch;
action_epoch = parameters.action_epoch;
sample_epoch = parameters.sample_epoch;
p_value_thresh = 0.05;

%% CNO
right_control_trial_idx = DREADD_activity.cno.trial_type.right_control_trial_idx;
left_control_trial_idx = DREADD_activity.cno.trial_type.left_control_trial_idx;
right_control_trial_choice = DREADD_activity.cno.trial_type.right_control_trial_choice;
left_control_trial_choice = DREADD_activity.cno.trial_type.left_control_trial_choice;

region_num = 1;
num_sample_cells{region_num} = nan(size(cno_all_resp,2),size(cno_all_resp{2},2));
num_choice_cells{region_num} = nan(size(cno_all_resp,2),size(cno_all_resp{2},2));
num_action_cells{region_num} = nan(size(cno_all_resp,2),size(cno_all_resp{2},2));

for animal_num = 1:size(cno_all_resp,2)
    
    for session_num = 1:size(cno_all_resp{animal_num},2)
        
        right_stim_activity = cno_all_resp{animal_num}{session_num}{region_num}(:,right_control_trial_idx{animal_num}{session_num},sample_epoch);
        left_stim_activity = cno_all_resp{animal_num}{session_num}{region_num}(:,left_control_trial_idx{animal_num}{session_num},sample_epoch);
        
        right_choice_activity = cno_all_resp{animal_num}{session_num}{region_num}(:,right_control_trial_choice{animal_num}{session_num},choice_epoch);
        left_choice_activity = cno_all_resp{animal_num}{session_num}{region_num}(:,left_control_trial_choice{animal_num}{session_num},choice_epoch);
        
        right_action_activity = cno_all_resp{animal_num}{session_num}{region_num}(:,right_control_trial_choice{animal_num}{session_num},action_epoch);
        left_action_activity = cno_all_resp{animal_num}{session_num}{region_num}(:,left_control_trial_choice{animal_num}{session_num},action_epoch);
        
        mean_right_stim_activity = squeeze(nanmean(right_stim_activity,3));
        mean_left_stim_activity= squeeze(nanmean(left_stim_activity,3));
        mean_right_choice_activity = squeeze(nanmean(right_choice_activity,3));
        mean_left_choice_activity = squeeze(nanmean(left_choice_activity,3));
        mean_right_action_activity = squeeze(nanmean(right_action_activity,3));
        mean_left_action_activity = squeeze(nanmean(left_action_activity,3));
        
        for cell_num = 1:size(mean_right_stim_activity,1)
            sample_p_value(cell_num) = ranksum(mean_right_stim_activity(cell_num,:),mean_left_stim_activity(cell_num,:),'tail','both');
            choice_p_value(cell_num) = ranksum(mean_right_choice_activity(cell_num,:),mean_left_choice_activity(cell_num,:),'tail','both');
            action_p_value(cell_num) = ranksum(mean_right_action_activity(cell_num,:),mean_left_action_activity(cell_num,:),'tail','both');
        end
        
        % Find cells with significant encoding activity
        sig_sample_cells{animal_num}{session_num}{region_num} = find( sample_p_value < p_value_thresh);
        sig_choice_cells{animal_num}{session_num}{region_num} = find( choice_p_value < p_value_thresh);
        sig_action_cells{animal_num}{session_num}{region_num} = find( action_p_value < p_value_thresh);
        
        clear  sample_p_value choice_p_value action_p_value
        
        for cell_num = 1:size(mean_right_stim_activity,1)
            sample_p_value(cell_num) = ranksum(mean_right_stim_activity(cell_num,:),mean_left_stim_activity(cell_num,:),'tail','right');
            choice_p_value(cell_num) = ranksum(mean_right_choice_activity(cell_num,:),mean_left_choice_activity(cell_num,:),'tail','right');
            action_p_value(cell_num) = ranksum(mean_right_action_activity(cell_num,:),mean_left_action_activity(cell_num,:),'tail','right');
        end
        
        % Find cells with significant encoding activity
        sig_sample_right_cells{animal_num}{session_num}{region_num} = find(sample_p_value < p_value_thresh);
        sig_choice_right_cells{animal_num}{session_num}{region_num} = find(choice_p_value < p_value_thresh);
        sig_action_right_cells{animal_num}{session_num}{region_num} = find(action_p_value < p_value_thresh);

        % Get fraction of cells with sig activity
        num_sample_cells{region_num}(animal_num,session_num) = length(sig_sample_cells{animal_num}{session_num}{region_num})./length(sample_p_value);
        num_choice_cells{region_num}(animal_num,session_num) = length(sig_choice_cells{animal_num}{session_num}{region_num})./length(choice_p_value);
        num_action_cells{region_num}(animal_num,session_num) = length(sig_action_cells{animal_num}{session_num}{region_num})./length(action_p_value);
        
        clear  sample_p_value choice_p_value action_p_value
    end
end

task_encoding.cno.sig_sample_cells = sig_sample_cells;
task_encoding.cno.sig_choice_cells = sig_choice_cells;
task_encoding.cno.sig_action_cells = sig_action_cells;

task_encoding.cno.sig_sample_right_cells = sig_sample_right_cells;
task_encoding.cno.sig_choice_right_cells = sig_choice_right_cells;
task_encoding.cno.sig_action_right_cells = sig_action_right_cells;

task_encoding.cno.num_sample_cells = num_sample_cells;
task_encoding.cno.num_choice_cells = num_choice_cells;
task_encoding.cno.num_action_cells = num_action_cells;

%% Saline
right_control_trial_idx = DREADD_activity.saline.trial_type.right_control_trial_idx;
left_control_trial_idx = DREADD_activity.saline.trial_type.left_control_trial_idx;
right_control_trial_choice = DREADD_activity.saline.trial_type.right_control_trial_choice;
left_control_trial_choice = DREADD_activity.saline.trial_type.left_control_trial_choice;

region_num = 1;
num_sample_cells{region_num} = nan(size(saline_all_resp,2),size(saline_all_resp{2},2));
num_choice_cells{region_num} = nan(size(saline_all_resp,2),size(saline_all_resp{2},2));
num_action_cells{region_num} = nan(size(saline_all_resp,2),size(saline_all_resp{2},2));
for animal_num = 1:size(saline_all_resp,2)
    for session_num = 1:size(saline_all_resp{animal_num},2)
        
        right_stim_activity = saline_all_resp{animal_num}{session_num}{region_num}(:,right_control_trial_idx{animal_num}{session_num},sample_epoch);
        left_stim_activity = saline_all_resp{animal_num}{session_num}{region_num}(:,left_control_trial_idx{animal_num}{session_num},sample_epoch);
        
        right_choice_activity = saline_all_resp{animal_num}{session_num}{region_num}(:,right_control_trial_choice{animal_num}{session_num},choice_epoch);
        left_choice_activity = saline_all_resp{animal_num}{session_num}{region_num}(:,left_control_trial_choice{animal_num}{session_num},choice_epoch);
        
        right_action_activity = saline_all_resp{animal_num}{session_num}{region_num}(:,right_control_trial_choice{animal_num}{session_num},action_epoch);
        left_action_activity = saline_all_resp{animal_num}{session_num}{region_num}(:,left_control_trial_choice{animal_num}{session_num},action_epoch);
        
        mean_right_stim_activity = squeeze(nanmean(right_stim_activity,3));
        mean_left_stim_activity= squeeze(nanmean(left_stim_activity,3));
        mean_right_choice_activity = squeeze(nanmean(right_choice_activity,3));
        mean_left_choice_activity = squeeze(nanmean(left_choice_activity,3));
        mean_right_action_activity = squeeze(nanmean(right_action_activity,3));
        mean_left_action_activity = squeeze(nanmean(left_action_activity,3));
        
        for cell_num = 1:size(mean_right_stim_activity,1)
            sample_p_value(cell_num) = ranksum(mean_right_stim_activity(cell_num,:),mean_left_stim_activity(cell_num,:),'tail','both');
            choice_p_value(cell_num) = ranksum(mean_right_choice_activity(cell_num,:),mean_left_choice_activity(cell_num,:),'tail','both');
            action_p_value(cell_num) = ranksum(mean_right_action_activity(cell_num,:),mean_left_action_activity(cell_num,:),'tail','both');
        end
        
        % Find cells with significant encoding activity
        sig_sample_cells{animal_num}{session_num}{region_num} = find( sample_p_value < p_value_thresh);
        sig_choice_cells{animal_num}{session_num}{region_num} = find( choice_p_value < p_value_thresh);
        sig_action_cells{animal_num}{session_num}{region_num} = find( action_p_value < p_value_thresh);
        
        clear  sample_p_value choice_p_value action_p_value
        
        for cell_num = 1:size(mean_right_stim_activity,1)
            sample_p_value(cell_num) = ranksum(mean_right_stim_activity(cell_num,:),mean_left_stim_activity(cell_num,:),'tail','right');
            choice_p_value(cell_num) = ranksum(mean_right_choice_activity(cell_num,:),mean_left_choice_activity(cell_num,:),'tail','right');
            action_p_value(cell_num) = ranksum(mean_right_action_activity(cell_num,:),mean_left_action_activity(cell_num,:),'tail','right');
        end
        
        % Find cells with significant encoding activity
        sig_sample_right_cells{animal_num}{session_num}{region_num} = find( sample_p_value < p_value_thresh);
        sig_choice_right_cells{animal_num}{session_num}{region_num} = find( choice_p_value < p_value_thresh);
        sig_action_right_cells{animal_num}{session_num}{region_num} = find( action_p_value < p_value_thresh);
        
        % Get fraction of cells with sig activity
        num_sample_cells{region_num}(animal_num,session_num) = length(sig_sample_cells{animal_num}{session_num}{region_num}) ./length(sample_p_value);
        num_choice_cells{region_num}(animal_num,session_num) = length(sig_choice_cells{animal_num}{session_num}{region_num})./length(choice_p_value);
        num_action_cells{region_num}(animal_num,session_num) = length(sig_action_cells{animal_num}{session_num}{region_num})./length(action_p_value);
        
        clear  sample_p_value choice_p_value action_p_value
    end
end

task_encoding.saline.sig_sample_cells = sig_sample_cells;
task_encoding.saline.sig_choice_cells = sig_choice_cells;
task_encoding.saline.sig_action_cells = sig_action_cells;
task_encoding.saline.sig_sample_right_cells = sig_sample_right_cells;
task_encoding.saline.sig_choice_right_cells = sig_choice_right_cells;
task_encoding.saline.sig_action_right_cells = sig_action_right_cells;
task_encoding.saline.num_sample_cells = num_sample_cells;
task_encoding.saline.num_choice_cells = num_choice_cells;
task_encoding.saline.num_action_cells = num_action_cells;
end