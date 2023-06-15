%% Plot decoding accuracy after getting AUC
function plot_auc_choice_stim_population(parameters,ROC_population_data,test_epoch)
concat_ROC_data = ROC_population_data.concat_auc;
concat_shuffled_auc_data = ROC_population_data.concat_shuffled_auc_data;
new_region_idx = parameters.new_region_idx;
fs_image = parameters.fs_image;
baseline_epoch = parameters.ITI_epoch;

%% Adjust by correcting against the shuffled baseline
baseline_epoch = parameters.ITI_epoch; %ITI
for session_num = 1:3
    for region_num = 1:8
        concat_ROC_data{session_num}{region_num} = concat_ROC_data{session_num}{region_num} - nanmean(concat_shuffled_auc_data{session_num}{region_num}(:,baseline_epoch),2);
        concat_shuffled_auc_data{session_num}{region_num} = concat_shuffled_auc_data{session_num}{region_num} - nanmean(concat_shuffled_auc_data{session_num}{region_num}(:,baseline_epoch),2);
    end
end

%% Test for significance between expert and naive at peak sensory diff
rng(2023)
for region_num = 1:8
    
    expert_data = concat_ROC_data{3}{region_num}(:,test_epoch);
    interm_data = concat_ROC_data{2}{region_num}(:,test_epoch);
    naive_data = concat_ROC_data{1}{region_num}(:,test_epoch);
    
    % Correlation bootstrap
    clear expert interm naive
    [~,max_idx(region_num)] = max( abs(nanmean(expert_data,1)-nanmean(naive_data,1)) );
    expert = expert_data(:,max_idx(region_num));
    interm = interm_data(:,max_idx(region_num));
    naive = naive_data(:,max_idx(region_num));
    for shuffle_num = 1:1000
        rand_idx_a = randsample(1:length(naive),length(naive),'true');
        rand_idx_b = randsample(1:length(interm),length(interm),'true');
        rand_idx_c = randsample(1:length(expert),length(expert),'true');
        
        [r_value(region_num,shuffle_num),~] = corr([1 2 3]',[nanmean(naive(rand_idx_a)) nanmean(interm(rand_idx_b)) nanmean(expert(rand_idx_c))]','Type','Pearson');
    end
end

% Bootstrap correlation
p_value = (sum(r_value < 0 ,2)./1000)';
p_value = p_value(parameters.new_region_idx);

% FDR Correction
clear fdr_p_value_decoding
for hide = 1
    % FDR correction
    q = 0.05;
    unroll = p_value;
    [a,b] = sort(unroll);
    [~,r] = sort(b);
    q_threshold = r/length(p_value) * q;
    fdr_p_value = unroll < q_threshold;
    q_threshold = r/length(p_value) * 0.01;
    fdr_p_value = fdr_p_value +(unroll < q_threshold);
    q_threshold = r/length(p_value) * 0.001;
    fdr_p_value = fdr_p_value +(unroll < q_threshold);
end

%% plot AUC trace
session_map(3,:) = [0.25 0.25 0.25];
session_map(2,:) = [0.5 0.5 0.5];
session_map(1,:) = [0.75 0.75 0.75];
figure('Position',[200,800,800,400],'Color','w','Color','white','DefaultAxesFontSize',14)
region_name = {'ALM','M1a','M1p','S1fl','vS1','M2','RSC','PPC'};
for i = 1:8
    region_num = (parameters.new_region_idx(i));
    subplot(2,4,i)
    hold on
    
    for session_num = [1 3]
        plot_data = concat_ROC_data{session_num}{region_num}';
        se = nanstd(plot_data,[],2) ./ sqrt(sum(~isnan(plot_data(1,:))));
        x = [1:size(plot_data,1)];
        x2 = [x fliplr(x)];
        curve1_1 = nanmean(plot_data,2)' + se';
        curve1_2 = nanmean(plot_data,2)' - se';
        in_between1 = [curve1_1 ,fliplr(curve1_2)];
        hold on
        h1 = fill(x2,in_between1,session_map(session_num,:),'LineStyle','none');
        set(h1,'facealpha',0.2)
        plot(nanmean(plot_data,2),'Color',session_map(session_num,:),'LineWidth',2)
        line([1.*fs_image,1.*fs_image],[-4,40],'Color',[0.25,0.25,0.25],'LineStyle',':')
        line([2.*fs_image,2.*fs_image],[-4,40],'Color',[0.25,0.25,0.25],'LineStyle',':')
        line([4.*fs_image,4.*fs_image],[-4,40],'Color',[0.25,0.25,0.25],'LineStyle',':')
        ax = gca;
        ax.XTickLabel = {''};
        ax.YTickLabel = {''};
        xlim([round(0*fs_image),round(4.*fs_image)])
        ylim([-0.1 1])
        if region_num == 1
            xlabel('Time to go cue')
            ylabel('Decoding Accuracy')
            yticks([0 0.5 1.0])
            yticklabels([0 0.5 1.0])
        elseif region_num ==2
            xticks([round(1.*fs_image) round(2.*fs_image)])
            xticklabels([-3 -2])
        else
            xticks([])
        end
        if fdr_p_value(i) ==3
            title([region_name{i},' ***'])
        elseif fdr_p_value(i) ==2
            title([region_name{i},' **'])
        elseif fdr_p_value(i) ==1
            title([region_name{i},' *'])
        else
            title(region_name{i})
        end
    end
end

%% Plot bar graph across region of difference.
% Test for significance between expert and naive at peak sensory diff
rng(2023)
for region_num = 1:8
    
    expert_data = concat_ROC_data{3}{region_num}(:,test_epoch);
    interm_data = concat_ROC_data{2}{region_num}(:,test_epoch);
    naive_data = concat_ROC_data{1}{region_num}(:,test_epoch);
    
    % Find the peak difference in expert-naive
    [~,max_idx(region_num)] = max( abs(nanmean(expert_data,1)-nanmean(naive_data,1)) );
    
    expert = expert_data(:,max_idx(region_num));
    interm= interm_data(:,max_idx(region_num));
    naive = naive_data(:,max_idx(region_num));
    
    all_expert{region_num} = expert;
    all_interm{region_num} = interm;
    all_naive{region_num} = naive;
    
    mean_all_expert(region_num) = nanmean(expert);
    mean_all_interm(region_num) = nanmean(interm);
    mean_all_naive(region_num) = nanmean(naive);
    
    sem_all_expert(region_num) = nanstd(expert) ./ sqrt(length(expert));
    sem_all_interm(region_num) = nanstd(interm) ./ sqrt(length(interm));
    sem_all_naive(region_num) = nanstd(naive) ./ sqrt(length(naive));
    
end

mean_AUC_naive = mean_all_naive(new_region_idx);
mean_AUC_intermediate = mean_all_interm(new_region_idx);
mean_AUC_expert = mean_all_expert(new_region_idx);

se_AUC_naive = sem_all_naive(new_region_idx);
se_AUC_intermediate = sem_all_interm(new_region_idx);
se_AUC_expert = sem_all_expert(new_region_idx);

figure('Position',[200,800,400,150],'Color','w')
hold on
bar(1,mean_AUC_naive(1),'FaceColor',session_map(1,:),'EdgeColor','None');
bar(2,mean_AUC_intermediate(1),'FaceColor',session_map(2,:),'EdgeColor','None');
bar(3,mean_AUC_expert(1),'FaceColor',session_map(3,:),'EdgeColor','None');
bar(5,mean_AUC_naive(2),'FaceColor',session_map(1,:),'EdgeColor','None');
bar(6,mean_AUC_intermediate(2),'FaceColor',session_map(2,:),'EdgeColor','None');
bar(7,mean_AUC_expert(2),'FaceColor',session_map(3,:),'EdgeColor','None');
bar(9,mean_AUC_naive(3),'FaceColor',session_map(1,:),'EdgeColor','None');
bar(10,mean_AUC_intermediate(3),'FaceColor',session_map(2,:),'EdgeColor','None');
bar(11,mean_AUC_expert(3),'FaceColor',session_map(3,:),'EdgeColor','None');
bar(13,mean_AUC_naive(4),'FaceColor',session_map(1,:),'EdgeColor','None');
bar(14,mean_AUC_intermediate(4),'FaceColor',session_map(2,:),'EdgeColor','None');
bar(15,mean_AUC_expert(4),'FaceColor',session_map(3,:),'EdgeColor','None');
bar(17,mean_AUC_naive(5),'FaceColor',session_map(1,:),'EdgeColor','None');
bar(18,mean_AUC_intermediate(5),'FaceColor',session_map(2,:),'EdgeColor','None');
bar(19,mean_AUC_expert(5),'FaceColor',session_map(3,:),'EdgeColor','None');
bar(21,mean_AUC_naive(6),'FaceColor',session_map(1,:),'EdgeColor','None');
bar(22,mean_AUC_intermediate(6),'FaceColor',session_map(2,:),'EdgeColor','None');
bar(23,mean_AUC_expert(6),'FaceColor',session_map(3,:),'EdgeColor','None');
bar(25,mean_AUC_naive(7),'FaceColor',session_map(1,:),'EdgeColor','None');
bar(26,mean_AUC_intermediate(7),'FaceColor',session_map(2,:),'EdgeColor','None');
bar(27,mean_AUC_expert(7),'FaceColor',session_map(3,:),'EdgeColor','None');
bar(29,mean_AUC_naive(8),'FaceColor',session_map(1,:),'EdgeColor','None');
bar(30,mean_AUC_intermediate(8),'FaceColor',session_map(2,:),'EdgeColor','None');
bar(31,mean_AUC_expert(8),'FaceColor',session_map(3,:),'EdgeColor','None');
line([1,1],[mean_AUC_naive(1) - se_AUC_naive(1),mean_AUC_naive(1) + se_AUC_naive(1)],'LineWidth',1,'Color',session_map(1,:))
line([2,2],[mean_AUC_intermediate(1) - se_AUC_intermediate(1),mean_AUC_intermediate(1) + se_AUC_intermediate(1)],'LineWidth',1,'Color',session_map(2,:))
line([3,3],[mean_AUC_expert(1) - se_AUC_expert(1),mean_AUC_expert(1) + se_AUC_expert(1)],'LineWidth',1,'Color',session_map(3,:))
line([5,5],[mean_AUC_naive(2) - se_AUC_naive(2),mean_AUC_naive(2) + se_AUC_naive(2)],'LineWidth',1,'Color',session_map(1,:))
line([6,6],[mean_AUC_intermediate(2) - se_AUC_intermediate(2),mean_AUC_intermediate(2) + se_AUC_intermediate(2)],'LineWidth',1,'Color',session_map(2,:))
line([7,7],[mean_AUC_expert(2) - se_AUC_expert(2),mean_AUC_expert(2) + se_AUC_expert(2)],'LineWidth',1,'Color',session_map(3,:))
line([9,9],[mean_AUC_naive(3) - se_AUC_naive(3),mean_AUC_naive(3) + se_AUC_naive(3)],'LineWidth',1,'Color',session_map(1,:))
line([10,10],[mean_AUC_intermediate(3) - se_AUC_intermediate(3),mean_AUC_intermediate(3) + se_AUC_intermediate(3)],'LineWidth',1,'Color',session_map(2,:))
line([11,11],[mean_AUC_expert(3) - se_AUC_expert(3),mean_AUC_expert(3) + se_AUC_expert(3)],'LineWidth',1,'Color',session_map(3,:))
line([13,13],[mean_AUC_naive(4) - se_AUC_naive(4),mean_AUC_naive(4) + se_AUC_naive(4)],'LineWidth',1,'Color',session_map(1,:))
line([14,14],[mean_AUC_intermediate(4) - se_AUC_intermediate(4),mean_AUC_intermediate(4) + se_AUC_intermediate(4)],'LineWidth',1,'Color',session_map(2,:))
line([15,15],[mean_AUC_expert(4) - se_AUC_expert(4),mean_AUC_expert(4) + se_AUC_expert(4)],'LineWidth',1,'Color',session_map(3,:))
line([17,17],[mean_AUC_naive(5) - se_AUC_naive(5),mean_AUC_naive(5) + se_AUC_naive(5)],'LineWidth',1,'Color',session_map(1,:))
line([18,18],[mean_AUC_intermediate(5) - se_AUC_intermediate(5),mean_AUC_intermediate(5) + se_AUC_intermediate(5)],'LineWidth',1,'Color',session_map(2,:))
line([19,19],[mean_AUC_expert(5) - se_AUC_expert(5),mean_AUC_expert(5) + se_AUC_expert(5)],'LineWidth',1,'Color',session_map(3,:))
line([21,21],[mean_AUC_naive(6) - se_AUC_naive(6),mean_AUC_naive(6) + se_AUC_naive(6)],'LineWidth',1,'Color',session_map(1,:))
line([22,22],[mean_AUC_intermediate(6) - se_AUC_intermediate(6),mean_AUC_intermediate(6) + se_AUC_intermediate(6)],'LineWidth',1,'Color',session_map(2,:))
line([23,23],[mean_AUC_expert(6) - se_AUC_expert(6),mean_AUC_expert(6) + se_AUC_expert(6)],'LineWidth',1,'Color',session_map(3,:))
line([25,25],[mean_AUC_naive(7) - se_AUC_naive(7),mean_AUC_naive(7) + se_AUC_naive(7)],'LineWidth',1,'Color',session_map(1,:))
line([26,26],[mean_AUC_intermediate(7) - se_AUC_intermediate(7),mean_AUC_intermediate(7) + se_AUC_intermediate(7)],'LineWidth',1,'Color',session_map(2,:))
line([27,27],[mean_AUC_expert(7) - se_AUC_expert(7),mean_AUC_expert(7) + se_AUC_expert(7)],'LineWidth',1,'Color',session_map(3,:))
line([29,29],[mean_AUC_naive(8) - se_AUC_naive(8),mean_AUC_naive(8) + se_AUC_naive(8)],'LineWidth',1,'Color',session_map(1,:))
line([30,30],[mean_AUC_intermediate(8) - se_AUC_intermediate(8),mean_AUC_intermediate(8) + se_AUC_intermediate(8)],'LineWidth',1,'Color',session_map(2,:))
line([31,31],[mean_AUC_expert(8) - se_AUC_expert(8),mean_AUC_expert(8) + se_AUC_expert(8)],'LineWidth',1,'Color',session_map(3,:))
xlim([-1,33])
ax = gca;
ax.FontSize = 14;
xticks([])
ylim([0 1])
yticks([0 0.5 1])
ylabel('Decoding accuracy')
for region_num = 1:8
    if fdr_p_value(region_num) == 3
        text( ((region_num*2 -1) *2),1,'***');
    elseif fdr_p_value(region_num) == 2
        text( ((region_num*2 -1) *2),1,'**');
    elseif fdr_p_value(region_num) == 1
        text( ((region_num*2 -1) *2),1,'*');
        
    end
end
end