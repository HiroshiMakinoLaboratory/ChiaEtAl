%% Get activity level
function average_activity = get_mean_activity_choice_encoding(parameters,DREADD_activity,task_encoding_cells)
% Initialise data
cno_all_resp = DREADD_activity.cno.activity.all_resp;
saline_all_resp = DREADD_activity.saline.activity.all_resp;
delay_epoch = [round(parameters.fs_image*3.5):round(parameters.fs_image*4)];
left_choice_trial = DREADD_activity.cno.trial_type.left_control_trial_choice;
right_choice_trial = DREADD_activity.cno.trial_type.right_control_trial_choice;

region_num = 1;
fs_image = parameters.fs_image;

%% CNO
choice_cell_idx = task_encoding_cells.cno.sig_choice_right_cells;
for animal_num = 1:size(cno_all_resp,2)
    for session_num = 1:size(cno_all_resp{animal_num},2)
        % Average across cells and trials
        left_choice_activity{animal_num}{session_num}{region_num} = squeeze(nanmean(mean(cno_all_resp{animal_num}{session_num}{region_num}(choice_cell_idx{animal_num}{session_num}{region_num},left_choice_trial{animal_num}{session_num},:),1),2));
        right_choice_activity{animal_num}{session_num}{region_num} = squeeze(nanmean(mean(cno_all_resp{animal_num}{session_num}{region_num}(choice_cell_idx{animal_num}{session_num}{region_num},right_choice_trial{animal_num}{session_num},:),1),2));
        
        % Average across epoch and trials
        left_choice_mean_activity{animal_num}{session_num}{region_num} = squeeze(nanmean(mean(cno_all_resp{animal_num}{session_num}{region_num}(choice_cell_idx{animal_num}{session_num}{region_num},left_choice_trial{animal_num}{session_num},delay_epoch),2),3));
        right_choice_mean_activity{animal_num}{session_num}{region_num} = squeeze(nanmean(mean(cno_all_resp{animal_num}{session_num}{region_num}(choice_cell_idx{animal_num}{session_num}{region_num},right_choice_trial{animal_num}{session_num},delay_epoch),2),3));
        
        % Mean/mean of all cells
        mean_all_left_choice_activity(animal_num,session_num) = nanmean(left_choice_mean_activity{animal_num}{session_num}{region_num});
        mean_all_right_choice_activity(animal_num,session_num) = nanmean(right_choice_mean_activity{animal_num}{session_num}{region_num});
    end
end

average_activity.cno.cells.left_choice_mean_activity = left_choice_mean_activity;
average_activity.cno.cells.right_choice_mean_activity = right_choice_mean_activity;

% Concat all animals together
counter = 1;
for animal_num = 1:size(saline_all_resp,2)
    for session_num = 1:size(saline_all_resp{animal_num},2)
        concat_left_choice_activity(counter,:) = left_choice_activity{animal_num}{session_num}{region_num}';
        concat_right_choice_activity(counter,:) = right_choice_activity{animal_num}{session_num}{region_num}';
        counter = counter+1;
    end
end

% remove the zeros
mean_all_left_choice_activity(mean_all_left_choice_activity==0) = nan;
mean_all_right_choice_activity(mean_all_right_choice_activity==0) = nan;
average_activity.cno.concat_left_choice_activity = concat_left_choice_activity;
average_activity.cno.concat_right_choice_activity = concat_right_choice_activity;
average_activity.cno.mean_all_left_choice_activity = mean_all_left_choice_activity;
average_activity.cno.mean_all_right_choice_activity = mean_all_right_choice_activity;

%% Saline
left_choice_trial = DREADD_activity.saline.trial_type.left_control_trial_choice;
right_choice_trial = DREADD_activity.saline.trial_type.right_control_trial_choice;
choice_cell_idx = task_encoding_cells.saline.sig_choice_right_cells;
for animal_num = 1:size(saline_all_resp,2)
    for session_num = 1:size(saline_all_resp{animal_num},2)
        
        % Average across cells and trials
        left_choice_activity{animal_num}{session_num}{region_num} = squeeze(nanmean(mean(saline_all_resp{animal_num}{session_num}{region_num}(choice_cell_idx{animal_num}{session_num}{region_num},left_choice_trial{animal_num}{session_num},:),1),2));
        right_choice_activity{animal_num}{session_num}{region_num} = squeeze(nanmean(mean(saline_all_resp{animal_num}{session_num}{region_num}(choice_cell_idx{animal_num}{session_num}{region_num},right_choice_trial{animal_num}{session_num},:),1),2));
        
        % Average across epoch and trials
        
        left_choice_mean_activity{animal_num}{session_num}{region_num} = squeeze(nanmean(mean(saline_all_resp{animal_num}{session_num}{region_num}(choice_cell_idx{animal_num}{session_num}{region_num},left_choice_trial{animal_num}{session_num},delay_epoch),2),3));
        right_choice_mean_activity{animal_num}{session_num}{region_num} = squeeze(nanmean(mean(saline_all_resp{animal_num}{session_num}{region_num}(choice_cell_idx{animal_num}{session_num}{region_num},right_choice_trial{animal_num}{session_num},delay_epoch),2),3));
        
        mean_all_left_choice_activity(animal_num,session_num) = nanmean(left_choice_mean_activity{animal_num}{session_num}{region_num});
        mean_all_right_choice_activity(animal_num,session_num) = nanmean(right_choice_mean_activity{animal_num}{session_num}{region_num});
    end
end
average_activity.saline.cells.left_choice_mean_activity = left_choice_mean_activity;
average_activity.saline.cells.right_choice_mean_activity = right_choice_mean_activity;

% Concat all animals together
counter = 1;
for animal_num = 1:size(saline_all_resp,2)
    for session_num = 1:size(saline_all_resp{animal_num},2)
        concat_left_choice_activity(counter,:) = left_choice_activity{animal_num}{session_num}{region_num}';
        concat_right_choice_activity(counter,:) = right_choice_activity{animal_num}{session_num}{region_num}';
        counter = counter+1;
    end
end

% remove the zeros
mean_all_left_choice_activity(mean_all_left_choice_activity==0) = nan;
mean_all_right_choice_activity(mean_all_right_choice_activity==0) = nan;
average_activity.saline.concat_left_choice_activity = concat_left_choice_activity;
average_activity.saline.concat_right_choice_activity = concat_right_choice_activity;
average_activity.saline.mean_all_left_choice_activity = mean_all_left_choice_activity;
average_activity.saline.mean_all_right_choice_activity = mean_all_right_choice_activity;

%% Plot Scatter plot
data1 = average_activity.cno.concat_left_choice_activity';
data2 = average_activity.cno.concat_right_choice_activity';
data3 = average_activity.saline.concat_left_choice_activity';
data4 = average_activity.saline.concat_right_choice_activity';

% Plot epoch average activity differences
epoch_analyse = [round(fs_image*3.5):round(fs_image*4)];
cno_data_diff = nanmean(data2(epoch_analyse,:) - data1(epoch_analyse,:));
saline_data_diff = nanmean(data4(epoch_analyse,:) - data3(epoch_analyse,:));

% Get significance using bootstrap
rng(2022)
plot_data = nanmean((data2(epoch_analyse,:) - data1(epoch_analyse,:)) - (data4(epoch_analyse,:) - data3(epoch_analyse,:)),1);
for shuffle_num = 1:1000
    shuffle_all(shuffle_num) = nanmean(plot_data(randsample(length(plot_data),length(plot_data),'true')));
end
p_value = sum(shuffle_all > 0)./1000;
average_activity.significance.p_value = p_value;

figure('Position',[200,100,150,150],'Color','white','DefaultAxesFontSize',14);
hold on;
scatter(cno_data_diff,saline_data_diff,'filled','MarkerFaceColor','k','MarkerEdgeColor','none')
line([0,0.4],[0,0.4],'Color',[0.25,0.25,0.25]);
xlim([0 0.4]);
ylim([0 0.4]);
cross_x = cno_data_diff;
cross_y = saline_data_diff;
SEM_a = nanstd(cross_x)./sqrt(length(cross_x));
SEM_b = nanstd(cross_y)./sqrt(length(cross_y));
line([nanmean(cross_x)-SEM_a nanmean(cross_x)+SEM_a], [nanmean(cross_y) nanmean(cross_y)],'Color',[0.85 0.325 0.098],'LineWidth',2)
line([nanmean(cross_x) nanmean(cross_x)], [nanmean(cross_y)-SEM_b nanmean(cross_y)+SEM_b],'Color',[0.85 0.325 0.098],'LineWidth',2)
if p_value < 0.001
    title('***')
elseif p_value < 0.01
    title('**')
elseif p_value < 0.05
    title('*')
end
xticks([0 0.4])
yticks([0 0.4])
xlabel('Choice selectivity (CNO)')
ylabel('Choice selectivity (Saline)')
axis square
end