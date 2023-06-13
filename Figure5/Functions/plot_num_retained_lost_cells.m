%% Plot fracton of retained and eliminated neurons
function [prop_retained,prop_lost] = plot_num_retained_lost_cells(retained_elim_data)
for animal_num = 1:13
    % Get Fraction of cells for each animal
    prop_retained(animal_num) = nanmean(retained_elim_data.num_cells.total_num_retained_naive_coupling{animal_num}./ (retained_elim_data.num_cells.total_num_retained_naive_coupling{animal_num}+retained_elim_data.num_cells.total_num_lost_naive_coupling{animal_num}) );
    prop_lost(animal_num) = nanmean(retained_elim_data.num_cells.total_num_lost_naive_coupling{animal_num}./ (retained_elim_data.num_cells.total_num_retained_naive_coupling{animal_num}+retained_elim_data.num_cells.total_num_lost_naive_coupling{animal_num}) );
    
end

cMap = colormap(lines(5));
figure('Position',[200,100,150,150],'Color','white','DefaultAxesFontSize',14);
hold on;
bar(1,mean(prop_retained),'FaceColor',cMap(2,:),'EdgeColor','None','BarWidth',0.7);
bar(2,mean(prop_lost),'FaceColor',cMap(1,:),'EdgeColor','None','BarWidth',0.7);
errorbar(1,mean(prop_retained),nanstd(prop_retained)/sqrt(length(prop_retained)),'CapSize',0,'Color',cMap(2,:));
errorbar(2,mean(prop_lost),nanstd(prop_lost)/sqrt(length(prop_lost)),'CapSize',0,'Color',cMap(1,:));
ylim([0 0.8]);
xticks([])
xlim([0 3])
yticks([0 0.4 0.8])
ylabel('Coupling fraction')
end