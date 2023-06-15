%% Get trial-by-trial correlation between regions
function coordination_data = get_coordination_data(parameters,data)
fs_image = parameters.fs_image;

for animal_num = 1:7
    for date_num = 1:6
        clearvars -except animal_list animal_num animal_ID date_list date_num trial_by_trial_corr_sensory_both trial_by_trial_corr_sensory_left trial_by_trial_corr_sensory_right trial_by_trial_corr_choice_both trial_by_trial_corr_choice_left trial_by_trial_corr_choice_right ...
            trial_by_trial_corr_sensory_delay_both trial_by_trial_corr_sensory_delay_left trial_by_trial_corr_sensory_delay_right trial_by_trial_corr_choice_delay_both trial_by_trial_corr_choice_delay_left trial_by_trial_corr_choice_delay_right result_base_folder...
            test_slope_sensory test_slope_choice change_sensory change_choice file_name data coordination_data fs_image
        
        % Load.
        sensory_mode_GS = data.all_sensory_mode_GS{animal_num}{date_num};
        choice_mode_GS = data.all_choice_mode_GS{animal_num}{date_num};
        SessionData = data.all_SessionData{animal_num}{date_num};
        
        % Based on trial type.
        left_sensory_trial = find(SessionData.TrialTypes == 1);
        right_sensory_trial = find(SessionData.TrialTypes == 2);
        
        % Based on choice.
        left_choice_trial = find(SessionData.TrialTypes == 1 & SessionData.Outcomes == 1 | SessionData.TrialTypes == 2 & SessionData.Outcomes == 0);
        right_choice_trial = find(SessionData.TrialTypes == 2 & SessionData.Outcomes == 1 | SessionData.TrialTypes == 1 & SessionData.Outcomes == 0);
        
        for region_num1 = 1:8
            for region_num2 = 1:8
                trial_by_trial_corr_sensory_both(animal_num,date_num,region_num1,region_num2) = corr(nanmean(sensory_mode_GS{region_num1}(:,round(1.*fs_image):round(2.*fs_image)),2),nanmean(sensory_mode_GS{region_num2}(:,round(1.*fs_image):round(2.*fs_image)),2),'Type','Pearson');
                trial_by_trial_corr_sensory_left(animal_num,date_num,region_num1,region_num2) = corr(nanmean(sensory_mode_GS{region_num1}(left_sensory_trial,round(1.*fs_image):round(2.*fs_image)),2),nanmean(sensory_mode_GS{region_num2}(left_sensory_trial,round(1.*fs_image):round(2.*fs_image)),2),'Type','Pearson');
                trial_by_trial_corr_sensory_right(animal_num,date_num,region_num1,region_num2) = corr(nanmean(sensory_mode_GS{region_num1}(right_sensory_trial,round(1.*fs_image):round(2.*fs_image)),2),nanmean(sensory_mode_GS{region_num2}(right_sensory_trial,round(1.*fs_image):round(2.*fs_image)),2),'Type','Pearson');
                trial_by_trial_corr_choice_both(animal_num,date_num,region_num1,region_num2) = corr(nanmean(choice_mode_GS{region_num1}(:,round(3.*fs_image):round(4.*fs_image)),2),nanmean(choice_mode_GS{region_num2}(:,round(3.*fs_image):round(4.*fs_image)),2),'Type','Pearson');
                if ~isempty(left_choice_trial) & length(left_choice_trial) > 5
                    trial_by_trial_corr_choice_left(animal_num,date_num,region_num1,region_num2) = corr(nanmean(choice_mode_GS{region_num1}(left_choice_trial,round(3.*fs_image):round(4.*fs_image)),2),nanmean(choice_mode_GS{region_num2}(left_choice_trial,round(3.*fs_image):round(4.*fs_image)),2),'Type','Pearson');
                else
                    trial_by_trial_corr_choice_left(animal_num,date_num,region_num1,region_num2) = nan;
                end
                if ~isempty(right_choice_trial) & length(right_choice_trial) > 5
                    trial_by_trial_corr_choice_right(animal_num,date_num,region_num1,region_num2) = corr(nanmean(choice_mode_GS{region_num1}(right_choice_trial,round(3.*fs_image):round(4.*fs_image)),2),nanmean(choice_mode_GS{region_num2}(right_choice_trial,round(3.*fs_image):round(4.*fs_image)),2),'Type','Pearson');
                else
                    trial_by_trial_corr_choice_right(animal_num,date_num,region_num1,region_num2) = nan;
                end
            end
        end
    end
end

for animal_num = 1:7
    trial_by_trial_corr_sensory_both_naive(animal_num,:,:,:) = squeeze(trial_by_trial_corr_sensory_both(animal_num,1:2,:,:));
    trial_by_trial_corr_sensory_both_intermediate(animal_num,:,:,:) = squeeze(trial_by_trial_corr_sensory_both(animal_num,3:4,:,:));
    trial_by_trial_corr_sensory_both_expert(animal_num,:,:,:) = squeeze(trial_by_trial_corr_sensory_both(animal_num,5:6,:,:));
    trial_by_trial_corr_sensory_left_naive(animal_num,:,:,:) = squeeze(trial_by_trial_corr_sensory_left(animal_num,1:2,:,:));
    trial_by_trial_corr_sensory_left_intermediate(animal_num,:,:,:) = squeeze(trial_by_trial_corr_sensory_left(animal_num,3:4,:,:));
    trial_by_trial_corr_sensory_left_expert(animal_num,:,:,:) = squeeze(trial_by_trial_corr_sensory_left(animal_num,5:6,:,:));
    trial_by_trial_corr_sensory_right_naive(animal_num,:,:,:) = squeeze(trial_by_trial_corr_sensory_right(animal_num,1:2,:,:));
    trial_by_trial_corr_sensory_right_intermediate(animal_num,:,:,:) = squeeze(trial_by_trial_corr_sensory_right(animal_num,3:4,:,:));
    trial_by_trial_corr_sensory_right_expert(animal_num,:,:,:) = squeeze(trial_by_trial_corr_sensory_right(animal_num,5:6,:,:));
    
    trial_by_trial_corr_choice_both_naive(animal_num,:,:,:) = squeeze(trial_by_trial_corr_choice_both(animal_num,1:2,:,:));
    trial_by_trial_corr_choice_both_intermediate(animal_num,:,:,:) = squeeze(trial_by_trial_corr_choice_both(animal_num,3:4,:,:));
    trial_by_trial_corr_choice_both_expert(animal_num,:,:,:) = squeeze(trial_by_trial_corr_choice_both(animal_num,5:6,:,:));
    trial_by_trial_corr_choice_left_naive(animal_num,:,:,:) = squeeze(trial_by_trial_corr_choice_left(animal_num,1:2,:,:));
    trial_by_trial_corr_choice_left_intermediate(animal_num,:,:,:) = squeeze(trial_by_trial_corr_choice_left(animal_num,3:4,:,:));
    trial_by_trial_corr_choice_left_expert(animal_num,:,:,:) = squeeze(trial_by_trial_corr_choice_left(animal_num,5:6,:,:));
    trial_by_trial_corr_choice_right_naive(animal_num,:,:,:) = squeeze(trial_by_trial_corr_choice_right(animal_num,1:2,:,:));
    trial_by_trial_corr_choice_right_intermediate(animal_num,:,:,:) = squeeze(trial_by_trial_corr_choice_right(animal_num,3:4,:,:));
    trial_by_trial_corr_choice_right_expert(animal_num,:,:,:) = squeeze(trial_by_trial_corr_choice_right(animal_num,5:6,:,:));
    
end

trial_by_trial_corr_sensory_both_naive_concat = [];
trial_by_trial_corr_sensory_both_intermediate_concat = [];
trial_by_trial_corr_sensory_both_expert_concat = [];
trial_by_trial_corr_sensory_left_naive_concat = [];
trial_by_trial_corr_sensory_left_intermediate_concat = [];
trial_by_trial_corr_sensory_left_expert_concat = [];
trial_by_trial_corr_sensory_right_naive_concat = [];
trial_by_trial_corr_sensory_right_intermediate_concat = [];
trial_by_trial_corr_sensory_right_expert_concat = [];

trial_by_trial_corr_choice_both_naive_concat = [];
trial_by_trial_corr_choice_both_intermediate_concat = [];
trial_by_trial_corr_choice_both_expert_concat = [];
trial_by_trial_corr_choice_left_naive_concat = [];
trial_by_trial_corr_choice_left_intermediate_concat = [];
trial_by_trial_corr_choice_left_expert_concat = [];
trial_by_trial_corr_choice_right_naive_concat = [];
trial_by_trial_corr_choice_right_intermediate_concat = [];
trial_by_trial_corr_choice_right_expert_concat = [];

for animal_num = 1:7
    trial_by_trial_corr_sensory_both_naive_concat = cat(1,trial_by_trial_corr_sensory_both_naive_concat,squeeze(trial_by_trial_corr_sensory_both_naive(animal_num,:,:,:)));
    trial_by_trial_corr_sensory_both_intermediate_concat = cat(1,trial_by_trial_corr_sensory_both_intermediate_concat,squeeze(trial_by_trial_corr_sensory_both_intermediate(animal_num,:,:,:)));
    trial_by_trial_corr_sensory_both_expert_concat = cat(1,trial_by_trial_corr_sensory_both_expert_concat,squeeze(trial_by_trial_corr_sensory_both_expert(animal_num,:,:,:)));
    trial_by_trial_corr_sensory_left_naive_concat = cat(1,trial_by_trial_corr_sensory_left_naive_concat,squeeze(trial_by_trial_corr_sensory_left_naive(animal_num,:,:,:)));
    trial_by_trial_corr_sensory_left_intermediate_concat = cat(1,trial_by_trial_corr_sensory_left_intermediate_concat,squeeze(trial_by_trial_corr_sensory_left_intermediate(animal_num,:,:,:)));
    trial_by_trial_corr_sensory_left_expert_concat = cat(1,trial_by_trial_corr_sensory_left_expert_concat,squeeze(trial_by_trial_corr_sensory_left_expert(animal_num,:,:,:)));
    trial_by_trial_corr_sensory_right_naive_concat = cat(1,trial_by_trial_corr_sensory_right_naive_concat,squeeze(trial_by_trial_corr_sensory_right_naive(animal_num,:,:,:)));
    trial_by_trial_corr_sensory_right_intermediate_concat = cat(1,trial_by_trial_corr_sensory_right_intermediate_concat,squeeze(trial_by_trial_corr_sensory_right_intermediate(animal_num,:,:,:)));
    trial_by_trial_corr_sensory_right_expert_concat = cat(1,trial_by_trial_corr_sensory_right_expert_concat,squeeze(trial_by_trial_corr_sensory_right_expert(animal_num,:,:,:)));
    
    trial_by_trial_corr_choice_both_naive_concat = cat(1,trial_by_trial_corr_choice_both_naive_concat,squeeze(trial_by_trial_corr_choice_both_naive(animal_num,:,:,:)));
    trial_by_trial_corr_choice_both_intermediate_concat = cat(1,trial_by_trial_corr_choice_both_intermediate_concat,squeeze(trial_by_trial_corr_choice_both_intermediate(animal_num,:,:,:)));
    trial_by_trial_corr_choice_both_expert_concat = cat(1,trial_by_trial_corr_choice_both_expert_concat,squeeze(trial_by_trial_corr_choice_both_expert(animal_num,:,:,:)));
    trial_by_trial_corr_choice_left_naive_concat = cat(1,trial_by_trial_corr_choice_left_naive_concat,squeeze(trial_by_trial_corr_choice_left_naive(animal_num,:,:,:)));
    trial_by_trial_corr_choice_left_intermediate_concat = cat(1,trial_by_trial_corr_choice_left_intermediate_concat,squeeze(trial_by_trial_corr_choice_left_intermediate(animal_num,:,:,:)));
    trial_by_trial_corr_choice_left_expert_concat = cat(1,trial_by_trial_corr_choice_left_expert_concat,squeeze(trial_by_trial_corr_choice_left_expert(animal_num,:,:,:)));
    trial_by_trial_corr_choice_right_naive_concat = cat(1,trial_by_trial_corr_choice_right_naive_concat,squeeze(trial_by_trial_corr_choice_right_naive(animal_num,:,:,:)));
    trial_by_trial_corr_choice_right_intermediate_concat = cat(1,trial_by_trial_corr_choice_right_intermediate_concat,squeeze(trial_by_trial_corr_choice_right_intermediate(animal_num,:,:,:)));
    trial_by_trial_corr_choice_right_expert_concat = cat(1,trial_by_trial_corr_choice_right_expert_concat,squeeze(trial_by_trial_corr_choice_right_expert(animal_num,:,:,:)));
end

coordination_data.trial_by_trial_corr_sensory_both_naive_concat = trial_by_trial_corr_sensory_both_naive_concat;
coordination_data.trial_by_trial_corr_sensory_both_intermediate_concat = trial_by_trial_corr_sensory_both_intermediate_concat;
coordination_data.trial_by_trial_corr_sensory_both_expert_concat = trial_by_trial_corr_sensory_both_expert_concat;
coordination_data.trial_by_trial_corr_sensory_left_naive_concat = trial_by_trial_corr_sensory_left_naive_concat;
coordination_data.trial_by_trial_corr_sensory_left_intermediate_concat = trial_by_trial_corr_sensory_left_intermediate_concat;
coordination_data.trial_by_trial_corr_sensory_left_expert_concat = trial_by_trial_corr_sensory_left_expert_concat;
coordination_data.trial_by_trial_corr_sensory_right_naive_concat = trial_by_trial_corr_sensory_right_naive_concat;
coordination_data.trial_by_trial_corr_sensory_right_intermediate_concat = trial_by_trial_corr_sensory_right_intermediate_concat;
coordination_data.trial_by_trial_corr_sensory_right_expert_concat = trial_by_trial_corr_sensory_right_expert_concat;

coordination_data.trial_by_trial_corr_choice_both_naive_concat = trial_by_trial_corr_choice_both_naive_concat;
coordination_data.trial_by_trial_corr_choice_both_intermediate_concat = trial_by_trial_corr_choice_both_intermediate_concat;
coordination_data.trial_by_trial_corr_choice_both_expert_concat = trial_by_trial_corr_choice_both_expert_concat;
coordination_data.trial_by_trial_corr_choice_left_naive_concat = trial_by_trial_corr_choice_left_naive_concat;
coordination_data.trial_by_trial_corr_choice_left_intermediate_concat = trial_by_trial_corr_choice_left_intermediate_concat;
coordination_data.trial_by_trial_corr_choice_left_expert_concat = trial_by_trial_corr_choice_left_expert_concat;
coordination_data.trial_by_trial_corr_choice_right_naive_concat = trial_by_trial_corr_choice_right_naive_concat;
coordination_data.trial_by_trial_corr_choice_right_intermediate_concat = trial_by_trial_corr_choice_right_intermediate_concat;
coordination_data.trial_by_trial_corr_choice_right_expert_concat = trial_by_trial_corr_choice_right_expert_concat;
end