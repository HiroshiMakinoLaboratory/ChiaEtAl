%% Get regional convergence index
function plot_conv_idx_heatmap(parameters,convergence_idx)
region = 8; % 8 cortical regions
conv_idx_heatmap_data = convergence_idx.matrix_mean_concat_regional_prop_other;
new_region_idx = parameters.new_region_idx;

figure('Position',[200,200,545,140],'Color','white','DefaultAxesFontSize',14);
c_axis = [0.05 0.35];
cMap = colormap(magma);
for session_num = 1:3
    subplot(1,3,session_num)
    imagesc(conv_idx_heatmap_data{session_num}(new_region_idx,new_region_idx));
    colormap(cMap);
    caxis(c_axis);
    xticks([]);yticks([])
    axis equal
    axis off
    if session_num == 1
        title('naive')
    elseif session_num == 2
        title('Intermediate')
    elseif session_num == 3
        title('Expert')
    end
    c= colorbar;
    c.Ticks = [min(c_axis) max(c_axis)];
    c.TickLabels = [min(c_axis) max(c_axis)];
end

%% Intra vs Inter regional
% Get significance using 3 way ANOVA
intra_only =  convergence_idx.raw.intra_only;
inter_only = convergence_idx.raw.inter_only;
for session_num = 1:3
    concat_all_intra{session_num} = []; region_label_intra{session_num} = []; concat_all_inter{session_num} = []; region_label_inter{session_num} = [];
    
    for region_num = 1:8
        concat_all_intra{session_num} = [concat_all_intra{session_num} intra_only{session_num}{region_num}];
        concat_all_inter{session_num} = [concat_all_inter{session_num} inter_only{session_num}{region_num}];
        
        region_label_intra{session_num} = [region_label_intra{session_num} repelem(region_num,length(intra_only{session_num}{region_num}))];
        region_label_inter{session_num} = [region_label_inter{session_num} repelem(region_num,length(inter_only{session_num}{region_num}))];
    end
end

% Concat labels to do 3-way ANOVA
concat_all = [concat_all_intra{1} concat_all_intra{2} concat_all_intra{3} concat_all_inter{1} concat_all_inter{2} concat_all_inter{3}];
concat_all_region = [region_label_intra{1} region_label_intra{2}  region_label_intra{3} region_label_inter{1} region_label_inter{2}  region_label_inter{3} ];
concat_all_session = [repelem(1,length(region_label_intra{1})) repelem(2,length(region_label_intra{2})) repelem(3,length(region_label_intra{3}))...
    repelem(1,length(region_label_inter{1})) repelem(2,length(region_label_inter{2})) repelem(3,length(region_label_inter{3}))];
concat_all_interintra = [repelem(1,length([concat_all_intra{1} concat_all_intra{2} concat_all_intra{3}])) repelem(12,length([concat_all_inter{1} concat_all_inter{2} concat_all_inter{3}]))];
anovan(concat_all,{concat_all_interintra concat_all_region concat_all_session});

% Plot
mean_intra_only = convergence_idx.inter_intra.mean_intra_only;
sem_intra_only = convergence_idx.inter_intra.sem_intra_only;
mean_inter_only = convergence_idx.inter_intra.mean_inter_only;
sem_inter_only = convergence_idx.inter_intra.sem_inter_only;

% Plot inter-regional conv. idx
figure('Position',[200,200,250,120],'Color','white','DefaultAxesFontSize',14);
y_lim = [0 0.65]; sessions_color_map = [222 222 222; 145 145 145 ; 52 52 52]./255;
for session_num = 1:3
    plot([1:8],mean_intra_only{session_num}(new_region_idx),'-','Color',sessions_color_map(session_num,:),'LineWidth',1,'MarkerFaceColor',sessions_color_map(session_num,:))
    hold on;
    box(gca,'off')
    xlim([0.5 region+0.5])
    xticks([1:8]);
    xticklabels([])
    ylim(y_lim)
    yticks([y_lim(1) y_lim(2)])
    errorbar([1:8],mean_intra_only{session_num}(new_region_idx),sem_intra_only{session_num},'CapSize',0,'Color',sessions_color_map(session_num,:),'LineWidth',2);
    hold on;
end

% Plot inter-regional conv. idx
for session_num = 1:3
    plot([1:8],mean_inter_only{session_num}(new_region_idx),':','Color',sessions_color_map(session_num,:),'LineWidth',1,'MarkerFaceColor',sessions_color_map(session_num,:))
    hold on;
    box(gca,'off')
    xlim([0.5 region+0.5])
    xticks([1:8]);
    xticklabels([])
    ylim(y_lim)
    yticks([y_lim(1) y_lim(2)])
    e = errorbar([1:8],mean_inter_only{session_num}(new_region_idx),sem_inter_only{session_num},'CapSize',0,'Color',sessions_color_map(session_num,:),'LineWidth',1);
    e.LineStyle = 'none';
    hold on;
    ylabel('mean convergence index')
end
end