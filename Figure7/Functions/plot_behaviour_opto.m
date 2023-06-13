%% Plot figure for optogenetic performance
function opto_behaviour = plot_behaviour_opto(ppc_correct_rate,vs1_correct_rate,session_idx)
%% Do bootstrap for significance

rng(2023)
a = squeeze(ppc_correct_rate(session_idx,1,:));
b = squeeze(ppc_correct_rate(session_idx,2,:));
for shuffle_num = 1:1000
    rand_idx = randsample(length(a(:)),length(a(:)),'true');
    shuffle_all(shuffle_num) = nanmean(a(rand_idx) - b(rand_idx));
end
p_value(1) = sum(shuffle_all<0)./1000;

a = squeeze(vs1_correct_rate(session_idx,1,:));
b = squeeze(vs1_correct_rate(session_idx,2,:));
for shuffle_num = 1:1000
    rand_idx = randsample(length(a(:)),length(a(:)),'true');
    shuffle_all(shuffle_num) = nanmean(a(rand_idx) - b(rand_idx));
end
p_value(2) = sum(shuffle_all<0)./1000;

%% Plot figures
% PPC
figure('Position',[200,200,150,150],'Color','white','DefaultAxesFontSize',14);
hold on
a = squeeze(ppc_correct_rate(session_idx,1,:));
b = squeeze(ppc_correct_rate(session_idx,2,:));
scatter(b(:),a(:),12,'filled','MarkerFaceColor','k','MarkerEdgeColor','k')
xlim([0.5 1.0])
ylim([0.5 1.0])
line([0,1],[0 1],'Color',[0.25,0.25,0.25]);
xticks([0.5 0.75 1.0])
yticks([0.5 0.75 1.0])
plot_cross_a = b(:);
plot_cross_b = a(:);
SEM_a = nanstd(plot_cross_a)./sqrt(length(plot_cross_a));
SEM_b = nanstd(plot_cross_b)./sqrt(length(plot_cross_b));
line([nanmean(plot_cross_a)-SEM_a nanmean(plot_cross_a)+SEM_a], [nanmean(plot_cross_b) nanmean(plot_cross_b)],'Color',[0.85 0.325 0.098],'LineWidth',2)
line([nanmean(plot_cross_a) nanmean(plot_cross_a)], [nanmean(plot_cross_b)-SEM_b nanmean(plot_cross_b)+SEM_b],'Color',[0.85 0.325 0.098],'LineWidth',2)
axis square
if p_value(1) < 0.001
    title('PPC ***')
elseif p_value(1) < 0.01
    title('PPC **')
elseif p_value(1) < 0.05
    title('PPC *')
else
    title('PPC')
end
xticklabels([50 75 100])
yticklabels([50 75 100])
ylabel('Correct (%) Control Trials')
xlabel('Correct (%) Opto. Trials')

% vS1
figure('Position',[200,200,150,150],'Color','white','DefaultAxesFontSize',14);
hold on
a = squeeze(vs1_correct_rate(session_idx,1,:));
b = squeeze(vs1_correct_rate(session_idx,2,:));
scatter(b(:),a(:),12,'filled','MarkerFaceColor','k','MarkerEdgeColor','k')
xlim([0.5 1.0])
ylim([0.5 1.0])
line([0,1],[0 1],'Color',[0.25,0.25,0.25]);
xticks([0.5 0.75 1.0])
yticks([0.5 0.75 1.0])
axis square
plot_cross_a = b(:);
plot_cross_b = a(:);
SEM_a = nanstd(plot_cross_a)./sqrt(length(plot_cross_a));
SEM_b = nanstd(plot_cross_b)./sqrt(length(plot_cross_b));
line([nanmean(plot_cross_a)-SEM_a nanmean(plot_cross_a)+SEM_a], [nanmean(plot_cross_b) nanmean(plot_cross_b)],'Color',[0.85 0.325 0.098],'LineWidth',2)
line([nanmean(plot_cross_a) nanmean(plot_cross_a)], [nanmean(plot_cross_b)-SEM_b nanmean(plot_cross_b)+SEM_b],'Color',[0.85 0.325 0.098],'LineWidth',2)
if p_value(2) < 0.001
    title('vS1 ***')
elseif p_value(2) < 0.01
    title('vS1 **')
elseif p_value(2) < 0.05
    title('vS1 *')
else 
    title('vS1')
end
xticklabels([50 75 100])
yticklabels([50 75 100])
ylabel('Correct (%) Control Trials')
xlabel('Correct (%) Opto. Trials')
opto_behaviour.significance.p_value = p_value;
end