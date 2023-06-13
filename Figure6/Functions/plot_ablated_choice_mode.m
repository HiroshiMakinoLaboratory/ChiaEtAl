%% Plot choice mode after selective ablation of regional functional coupling
function plot_ablated_choice_mode(parameters,ablated_data)
% Compare original vs control
fs_image = parameters.fs_image;

% Initialise parameters
session_num = 1;
cMap = colormap([ [0 0.4470 0.7410] ;[0.5 0.5 0.5] ]);
x_lim = [round(fs_image*2),round(fs_image*4)]; y_lim = [-2 2];
region_idx = [1 6 7 8]; % Regions ALM, vS1, RSC and PPC
without_resp_epoch = parameters.pre_resp_epoch;
stim_offset = parameters.sample_offset;

concat_original_right_choice_mode = ablated_data.choice_mode.concat_original_right_choice_mode;
concat_original_left_choice_mode = ablated_data.choice_mode.concat_original_left_choice_mode;
concat_control_right_choice_mode = ablated_data.choice_mode.concat_control_right_choice_mode;
concat_control_left_choice_mode= ablated_data.choice_mode.concat_control_left_choice_mode;

counter = 1;
figure('Position',[200,100,730,700],'Color','white','DefaultAxesFontSize',14);
for region_num = region_idx
    for other_region_num = region_idx
        
        concat_diff_data_1=[];concat_diff_data_2=[];
        
        clear mean_diff_data
        valid_animal_idx = find(abs(cellfun(@isempty,concat_original_right_choice_mode{session_num}{region_num,other_region_num}{1})-1));
        for animal_num = valid_animal_idx
            % Process each shuffle independently and concat all together
            for shuffle_num = 1:size(concat_original_right_choice_mode{session_num}{region_num,other_region_num}{1}{animal_num},3)
                
                original_data_right = concat_original_right_choice_mode{session_num}{region_num,other_region_num}{1}{animal_num}(:,without_resp_epoch,shuffle_num);
                original_data_left = concat_original_left_choice_mode{session_num}{region_num,other_region_num}{1}{animal_num}(:,without_resp_epoch,shuffle_num);
                
                ablated_data_right = concat_control_right_choice_mode{session_num}{region_num,other_region_num}{1}{animal_num}(:,without_resp_epoch,shuffle_num);
                ablated_data_left = concat_control_left_choice_mode{session_num}{region_num,other_region_num}{1}{animal_num}(:,without_resp_epoch,shuffle_num);
                
                % Average across trials
                original_data_right = nanmean(original_data_right,1);
                original_data_left = nanmean(original_data_left,1);
                
                ablated_data_right = nanmean(ablated_data_right,1);
                ablated_data_left = nanmean(ablated_data_left,1);
                
                % Baseline of each data
                offset_original_right  = nanmean(original_data_right(:,stim_offset),2);
                offset_original_left  = nanmean(original_data_left(:,stim_offset),2);
                
                offset_ablated_right  = nanmean(ablated_data_right(:,stim_offset),2);
                offset_ablated_left = nanmean(ablated_data_left(:,stim_offset),2);
                
                % Subtract baseline
                original_data_right = original_data_right - offset_original_right;
                original_data_left = original_data_left - offset_original_left;
                
                ablated_data_right = ablated_data_right - offset_ablated_right;
                ablated_data_left = ablated_data_left - offset_ablated_left;
                
                % Normalise by std so that it will be comparable across regions
                norm_original_data_right = original_data_right ./nanstd([original_data_right ablated_data_right]);
                norm_original_data_left = original_data_left ./nanstd([original_data_right ablated_data_right]);
                
                norm_ablated_data_right = ablated_data_right ./nanstd([original_data_right ablated_data_right]);
                norm_ablated_data_left = ablated_data_left ./nanstd([original_data_right ablated_data_right ]);
                
                % Get the difference
                diff_data_1(:,shuffle_num) = ((norm_original_data_right ))';
                diff_data_2(:,shuffle_num) = ((norm_ablated_data_right ))';
            end
            mean_diff_data_1 = nanmean(diff_data_1,2);
            mean_diff_data_2 = nanmean(diff_data_2,2);
            
            concat_diff_data_1 = [concat_diff_data_1 mean_diff_data_1];
            concat_diff_data_2 = [concat_diff_data_2 mean_diff_data_2];
            
            clear diff_data_1 diff_data_2
        end
        
        % Average across animals/sessions
        concat_mean_diff_data_1 = squeeze(nanmean(concat_diff_data_1,2));
        concat_mean_diff_data_2 = squeeze(nanmean(concat_diff_data_2,2));
        
        subplot(length(region_idx),length(region_idx),counter)
        
        % Control
        x = [1:size(concat_diff_data_1,1)];
        x2 = [x fliplr(x)];
        se = nanstd(concat_diff_data_1,[],2) ./ sqrt(sum(~isnan(concat_diff_data_1(1,:))));
        curve1_1 = nanmean(concat_diff_data_1,2)' + se';
        curve1_2 = nanmean(concat_diff_data_1,2)' - se';
        in_between1 = [curve1_1 ,fliplr(curve1_2)];
        hold on
        h1 = fill(x2,in_between1,cMap(2,:),'LineStyle','none','LineWidth',1);
        set(h1,'facealpha',0.2)
        hold on;
        plot(concat_mean_diff_data_1','Color',cMap(2,:));
        
        % Ablated
        x = [1:size(concat_diff_data_2,1)];
        x2 = [x fliplr(x)];
        se = nanstd(concat_diff_data_2,[],2) ./ sqrt(sum(~isnan(concat_diff_data_2(1,:))));
        curve1_1 = nanmean(concat_diff_data_2,2)' + se';
        curve1_2 = nanmean(concat_diff_data_2,2)' - se';
        in_between1 = [curve1_1 ,fliplr(curve1_2)];
        hold on
        h1 = fill(x2,in_between1,cMap(1,:),'LineStyle','none');
        set(h1,'facealpha',0.2)
        hold on;
        plot(concat_mean_diff_data_2','Color',cMap(1,:),'LineWidth',1);
        if region_num == 1 & other_region_num == 1
            ylabel('Choice activity')
            xlabel('Time')
            xticks([0 2])
            yticks([-2 0 2])
        else
            xticks([])
            yticks([])
        end
        box(gca,'off')
        xlim([x_lim])
        ylim([y_lim])
        
        counter = counter +1;
    end
end
end