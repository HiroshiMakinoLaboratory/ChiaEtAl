%% Plot overtrained dataset
function plot_diff_overtrained(ppc_correct_rate)

hold on;
session_idx = [1:3];
overtrained_session_idx = [4:6];

%% Perform bootstrap for significance
a = squeeze(ppc_correct_rate(session_idx,1,:));
b = squeeze(ppc_correct_rate(overtrained_session_idx,1,:));
a = a(~isnan(a));
b = b(~isnan(b));
% Comparing trained and overtrained sessions
for shuffle_num = 1:1000
    shuffle_all(shuffle_num) = nanmean(randsample(a,length(a),'true')) - nanmean(randsample(b,length(b),'true'));
end
p_value = sum(shuffle_all > 0) ./ 1000;

% Plot figure
figure('Position',[200,800,150,150],'Color','w');
hold on
bar([nanmean(a),nanmean(b)],'EdgeColor','none','FaceColor','k','BarWidth',0.7)
e = errorbar([1 2],[nanmean(a) nanmean(b)  ],[nanstd(a)./sqrt(length(a)) nanstd(b)./sqrt(length(b))],'k','CapSize',0);
e.LineStyle = 'none';
xticks([])
ylim([0.5 1])
yticks([0.5 0.75 1.0])
yticklabels([50 75 100])
box(gca,'off')
if p_value(1) < 0.001
    title('PPC ***')
elseif p_value(1) < 0.01
    title('PPC **')
elseif p_value(1) < 0.05
    title('PPC *')
else
    title('')
end
ylabel('Performance (%) control trials')
xticks([0.8 2.2])
xticklabels({'trained','overtrained'})

%% Perform bootstrap for significance
a = squeeze(ppc_correct_rate(session_idx,2,:));
b = squeeze(ppc_correct_rate(overtrained_session_idx,2,:));
a = a(~isnan(a));
b = b(~isnan(b));
for shuffle_num = 1:1000
    shuffle_all(shuffle_num) = nanmean(randsample(a,length(a),'true')) - nanmean(randsample(b,length(b),'true'));
end
p_value = sum(shuffle_all > 0) ./ 1000;

% Plot figure
figure('Position',[200,800,150,150],'Color','w');
hold on
bar([nanmean(a),nanmean(b)],'EdgeColor','none','FaceColor','k','BarWidth',0.7)
e = errorbar([1 2],[nanmean(a) nanmean(b)  ],[nanstd(a)./sqrt(length(a)) nanstd(b)./sqrt(length(b))],'k','CapSize',0);
e.LineStyle = 'none';
xticks([])
ylim([0.5 1])
yticks([0.5 0.75 1.0])
yticklabels([50 75 100])
box(gca,'off')
if p_value(1) < 0.001
    title('PPC ***')
elseif p_value(1) < 0.01
    title('PPC **')
elseif p_value(1) < 0.05
    title('PPC *')
else
    title('')
end

ylabel('Performance (%) opto. trials')
xticks([0.8 2.2])
xticklabels({'trained','overtrained'})
end