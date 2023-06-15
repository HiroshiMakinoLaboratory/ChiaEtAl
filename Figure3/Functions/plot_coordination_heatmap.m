%% Plot heatmap for coordination
function plot_coordination_heatmap(parameters, data, sensory_or_choice, both_or_same, session_1, session_2)
new_region_idx = parameters.new_region_idx;

if strcmp(both_or_same,'both')
    data1 = data.(['trial_by_trial_corr_',sensory_or_choice,'_',both_or_same,'_',session_1,'_concat']);
    data2 = data.(['trial_by_trial_corr_',sensory_or_choice,'_',both_or_same,'_',session_2,'_concat']);
    
    figure('Position',[200,800,400,200],'Color','w')
    subplot(1,2,1)
    imagesc(squeeze(nanmean(data1(:,new_region_idx,new_region_idx))),[0,1])
    axis square
    ax = gca;
    ax.XTickLabel = {''};
    ax.YTickLabel = {''};
    colormap(magma)
    axis off
    colorbar
    squeeze(sum(~isnan(data1(:,new_region_idx,new_region_idx)),1))
    
    subplot(1,2,2)
    imagesc(squeeze(nanmean(data2(:,new_region_idx,new_region_idx))),[0,1])
    axis square
    ax = gca;
    ax.XTickLabel = {''};
    ax.YTickLabel = {''};
    colormap(magma)
    axis off
    colorbar
    squeeze(sum(~isnan(data2(:,new_region_idx,new_region_idx)),1))
    
elseif strcmp(both_or_same,'same')
    
    data1 = data.(['trial_by_trial_corr_',sensory_or_choice,'_left_',session_1,'_concat']);
    data2 = data.(['trial_by_trial_corr_',sensory_or_choice,'_right_',session_1,'_concat']);
    data3 = data.(['trial_by_trial_corr_',sensory_or_choice,'_left_',session_2,'_concat']);
    data4 = data.(['trial_by_trial_corr_',sensory_or_choice,'_right_',session_2,'_concat']);
    
    % Left + right choice
    figure('Position',[200,800,400,200],'Color','w')
    subplot(1,2,1)
    data_plot = [data1(:,[1,2,3,5,6,4,7,8],[1,2,3,5,6,4,7,8]);data2(:,[1,2,3,5,6,4,7,8],[1,2,3,5,6,4,7,8])];
    imagesc(squeeze(nanmean(data_plot)),[0,1])
    axis square
    ax = gca;
    ax.XTickLabel = {''};
    ax.YTickLabel = {''};
    colormap(magma)
    axis off
    colorbar
    
    subplot(1,2,2)
    data_plot = [data3(:,[1,2,3,5,6,4,7,8],[1,2,3,5,6,4,7,8]);data4(:,[1,2,3,5,6,4,7,8],[1,2,3,5,6,4,7,8])];
    imagesc(squeeze(nanmean(data_plot)),[0,1])
    axis square
    ax = gca;
    ax.XTickLabel = {''};
    ax.YTickLabel = {''};
    colormap(magma)
    axis off
    colorbar
end
end