%% Plot DREADD behaviour data
function [concat_data,p_value] = plot_behaviour_DREADD(CNO_correct_rate,Saline_correct_rate)
% Compare control trials
all_cno_control_trials = CNO_correct_rate(:,1,:);
all_cno_control_trials = all_cno_control_trials(:);
all_cno_control_trials = all_cno_control_trials(all_cno_control_trials~=0);

all_saline_control_trials = Saline_correct_rate(:,1,:);
all_saline_control_trials = all_saline_control_trials(:);
all_saline_control_trials = all_saline_control_trials(all_cno_control_trials~=0);

concat_data = [all_cno_control_trials all_saline_control_trials];

%% Bootstrap for significance
rng(2020)
for shuffle_num = 1:1000
    rand_idx = randsample(length(all_cno_control_trials),length(all_cno_control_trials),'true');
    shuffle_all(shuffle_num) = nanmean(all_cno_control_trials(rand_idx) - all_saline_control_trials(rand_idx));
end
p_value = sum(shuffle_all>0)/1000;

figure('Position',[200,200,170,140],'Color','white','DefaultAxesFontSize',14);
hold on
scatter(all_cno_control_trials,all_saline_control_trials,'filled','MarkerFaceColor','k','MarkerEdgeColor','k')
xlim([0.5 1.0])
ylim([0.5 1.0])
line([0,1],[0 1],'Color',[0.25,0.25,0.25]);
cross_x = all_cno_control_trials;
cross_y = all_saline_control_trials;
SEM_a = nanstd(cross_x)./sqrt(length(cross_x));
SEM_b = nanstd(cross_y)./sqrt(length(cross_y));
line([nanmean(cross_x)-SEM_a nanmean(cross_x)+SEM_a], [nanmean(cross_y) nanmean(cross_y)],'Color',[0.85 0.325 0.098],'LineWidth',2)
line([nanmean(cross_x) nanmean(cross_x)], [nanmean(cross_y)-SEM_b nanmean(cross_y)+SEM_b],'Color',[0.85 0.325 0.098],'LineWidth',2)
xticks([0.5 0.75 1.0])
yticks([0.5 0.75 1.0])
xlabel('Correct rate (CNO)')
ylabel('Correct rate (Saline)')
axis square
if p_value < 0.001
    title('***')
elseif p_value < 0.01
    title('**')
elseif p_value < 0.05
    title('*')
end
end