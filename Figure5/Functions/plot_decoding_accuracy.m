%% Plot regional decoding accuracy for retained vs lost
function plot_decoding_accuracy(parameters,decoding_accuracy)
new_region_idx = parameters.new_region_idx;

% Load data
region_concat_retained_auc = decoding_accuracy.region_concat_retained_auc;
region_concat_lost_auc = decoding_accuracy.region_concat_lost_auc;

figure('Position',[200,100,300,150],'Color','white','DefaultAxesFontSize',14);
for region_num = new_region_idx
    hold on;
    a = (abs(region_concat_retained_auc{2}{region_num}(~isnan(region_concat_retained_auc{2}{region_num})) -0.5) *2);
    b = (abs(region_concat_lost_auc{2}{region_num}(~isnan(region_concat_lost_auc{2}{region_num})) -0.5) *2);
    diff = a-b;
    
    bar(region_num,mean(diff),'FaceColor',[0 0 0],'EdgeColor','None','BarWidth',0.6);
    e = errorbar(region_num,mean(a-b),nanstd(a-b)/sqrt(length(a-b)),'CapSize',0,'Color',[0 0 0]);
    xlim([0 9]);
    ylim([-0.01 0.11])
    xticks([]);
    ylabel('Diff. Decoding accuracy (retained - eliminated)')
    xlabel('Regions in "post"')
    
    total_num_cells(region_num) = length(b);
end
end