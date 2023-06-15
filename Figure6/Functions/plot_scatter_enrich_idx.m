%% Plot animal-averaged enrichment index (Figure 5b)
function plot_scatter_enrich_idx(data1)
valid_idx = find((data1(:,1) >0 ));
naive_data = data1(valid_idx,1);
interm_data = data1(valid_idx,2);
expert_data = data1(valid_idx,3);

% Test significance using bootstrap
rng(2020)
for shuffle_num = 1:1000
    rand_idx = randsample(length(valid_idx),length(valid_idx),'true');
    shuffle_mean(shuffle_num) = nanmean(expert_data(rand_idx) - naive_data(rand_idx));
end
p_value = sum(shuffle_mean<0) ./ 1000;

% Plot
figure('Position',[200,200,200,170],'Color','white','DefaultAxesFontSize',14);
hold on
scatter(naive_data,expert_data,'filled','MarkerFaceColor','k','MarkerEdgeColor','k')
xlim([0.2 0.9])
ylim([0.2 0.9])
line([0,1],[0 1],'Color',[0.25,0.25,0.25]);
cross_x = naive_data;
cross_y = expert_data;
SEM_a = nanstd(cross_x)./sqrt(length(cross_x));
SEM_b = nanstd(cross_y)./sqrt(length(cross_y));
line([nanmean(cross_x)-SEM_a nanmean(cross_x)+SEM_a], [nanmean(cross_y) nanmean(cross_y)],'Color',[0.85 0.325 0.098],'LineWidth',2)
line([nanmean(cross_x) nanmean(cross_x)], [nanmean(cross_y)-SEM_b nanmean(cross_y)+SEM_b],'Color',[0.85 0.325 0.098],'LineWidth',2)
if p_value < 0.001
    title('***')
elseif p_value < 0.01
    title('**')
elseif p_value < 0.05
    title('*')
end
ylabel('Enrich. index (expert)')
xlabel('Enrich. index (naive)')
axis square
end