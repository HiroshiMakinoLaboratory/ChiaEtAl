%% Get global convergence index
function plot_global_conv_idx(convergence_idx)
plot_data(:,1)  = convergence_idx.global.concat_conv_idx{1};
plot_data(:,2)  = convergence_idx.global.concat_conv_idx{2};
plot_data(:,3)  = convergence_idx.global.concat_conv_idx{3};

%% Get significance using one-way ANOVA
p_value = anova1(plot_data,[],'off');

%% Plot
figure('Position',[200,1000,150,140],'Color','white','DefaultAxesFontSize',14);
hold on;
bar(1, nanmean(plot_data(:,1)),'FaceColor',[0.75 0.75 0.75],'EdgeColor','None','FaceAlpha',1,'BarWidth',0.7) 
bar(2, nanmean(plot_data(:,2)),'FaceColor',[0.5 0.5 0.5],'EdgeColor','None','FaceAlpha',1,'BarWidth',0.7) 
bar(3, nanmean(plot_data(:,3)),'FaceColor',[0.25 0.25 0.25],'EdgeColor','None','FaceAlpha',1,'BarWidth',0.7) 

errorbar(1, nanmean(plot_data(:,1)), (nanstd(plot_data(:,1))./sqrt(length(plot_data(:,1))) ) ,'LineWidth',2,'CapSize',0,'LineStyle','none','Color',[0.75 0.75 0.75] )
errorbar(2, nanmean(plot_data(:,2)), (nanstd(plot_data(:,2))./sqrt(length(plot_data(:,2))) ) ,'LineWidth',2,'CapSize',0,'LineStyle','none','Color',[0.5 0.5 0.5] )
errorbar(3, nanmean(plot_data(:,3)), (nanstd(plot_data(:,3))./sqrt(length(plot_data(:,3))) ) ,'LineWidth',2,'CapSize',0,'LineStyle','none','Color',[0.25 0.25 0.25] )
box(gca,'off');

ylabel('Convergence index')
xlim([0 4])
ylim([0 0.4])
end