%% Plot behavior curve
function plot_behavior_curve(behavior_performance)

plot_data = behavior_performance;
figure('Position',[200,200,170,150],'Color','white','DefaultAxesFontSize',14);
hold on

% Plot individual mouse behavior
for i = 1:7
    plot([1:20],plot_data(i,:),'Color',[0.8 0.8 0.8])
end

% Plot shaded area
plot_data = plot_data'
se = nanstd(plot_data,[],2) ./ sqrt(sum(~isnan(plot_data(1,:))));
x = [1:size(plot_data,1)];
x2 = [x fliplr(x)];
curve1_1 = nanmean(plot_data,2)' + se';
curve1_2 = nanmean(plot_data,2)' - se';
in_between1 = [curve1_1 ,fliplr(curve1_2)];
hold on
h1 = fill(x2,in_between1,[0.5 0.5 0.5],'LineStyle','none');
set(h1,'facealpha',0.4)

% Plot mean
plot(nanmean(plot_data'),'k','LineWidth',1.5)

ylim([0 100])
yticks([0:20:100])
ylabel('Correct (%)')
xlim([1 20])
xticks([0:5:20])
xlabel('Session')

end
