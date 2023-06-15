%% Analyse changes in GLM weights over learning
function GLM_coeff = get_GLM_coeff(parameters,data,stim_choice_delay)
all_GLM = data.all_GLM;
all_valid_cells = data.cell_idx.all_valid_cells;
all_region_labels = data.cell_idx.all_region_labels;
% Get predictor index for GlM coefficients

if strcmp(stim_choice_delay,'stimulus')
    left_predictor_idx = [1 2 3];
    right_predictor_idx = [4 5 6];
elseif strcmp(stim_choice_delay,'delay')
    left_predictor_idx = [7:11];
    right_predictor_idx = [12:16];
elseif strcmp(stim_choice_delay,'action')
    left_predictor_idx = [17:21];
    right_predictor_idx = [22:26];
end

for animal_num = 1:13
    for session_num = 1:3
        % Stimulus
        left_coeff = all_GLM{animal_num}{session_num}.coeff(:,left_predictor_idx);
        right_coeff = all_GLM{animal_num}{session_num}.coeff(:,right_predictor_idx);
        
        % Remove cells that are not significantly predicted by the GLM
        valid_cell_idx = all_valid_cells{animal_num}{session_num};
        left_coeff(setdiff(1:size(left_coeff,1),valid_cell_idx),:) = nan;
        right_coeff(setdiff(1:size(right_coeff,1),valid_cell_idx),:) = nan;
        
        % Take the strongest predictor
        max_left_coeff{animal_num}(:,session_num) = max(left_coeff,[],2);
        max_stim_coeff{animal_num}(:,session_num) = max(right_coeff,[],2);
        
        % Get average coeff
        animal_left_coeff(animal_num,session_num) = nanmean(max_left_coeff{animal_num}(:,session_num));
        animal_right_coeff(animal_num,session_num) = nanmean(max_stim_coeff{animal_num}(:,session_num));
        
        for region_num = 1:8
            region_idx = intersect(find(all_region_labels{animal_num}{session_num} == region_num) ,valid_cell_idx);
            if length(region_idx) > 0
                regional_animal_left_coeff{region_num}(animal_num,session_num) =  nanmean(max_left_coeff{animal_num}(region_idx,session_num));
                regional_animal_right_coeff{region_num}(animal_num,session_num) =  nanmean(max_stim_coeff{animal_num}(region_idx,session_num));
            else
                regional_animal_left_coeff{region_num}(animal_num,session_num) =  nan;
                regional_animal_right_coeff{region_num}(animal_num,session_num) =  nan;
            end
        end
    end
end

% Save variables
GLM_coeff.stim.animal_left_stim_coeff = animal_left_coeff;
GLM_coeff.stim.animal_right_stim_coeff = animal_right_coeff;
GLM_coeff.stim.regional_animal_left_stim_coeff = regional_animal_left_coeff;
GLM_coeff.stim.regional_animal_right_stim_coeff = regional_animal_right_coeff;

%% Get P value using Repeated measures ANOVA
a1 = animal_left_coeff(:,1);
a2 = animal_left_coeff(:,2);
a3 = animal_left_coeff(:,3);
t = table(a1,a2,a3,'VariableNames',{'a1','a2','a3'});
time = [1 2 3]';
rm = fitrm(t,'a1-a3 ~ 1','WithinDesign',time);
ranovatbl = ranova(rm);
temp = table2cell(ranovatbl(1,5));
p_value(1) = temp{1};

a1 = animal_right_coeff(:,1);
a2 = animal_right_coeff(:,2);
a3 = animal_right_coeff(:,3);
t = table(a1,a2,a3,'VariableNames',{'a1','a2','a3'});
time = [1 2 3]';
rm = fitrm(t,'a1-a3 ~ 1','WithinDesign',time);
ranovatbl = ranova(rm);
temp = table2cell(ranovatbl(1,5));
p_value(2) = temp{1};

GLM_coeff.significance.p_value = p_value;

%% Plot global fraction
% Plot left
figure('Position',[200,800,300,150],'Color','w')
subplot(1,2,1)
plot_data = animal_left_coeff;
bar([1],nanmean(plot_data(:,1)),'FaceColor',[0.75,0.75,0.75],'EdgeColor','None','FaceAlpha',1); hold on
bar([2],nanmean(plot_data(:,2)),'FaceColor',[0.5,0.5,0.5],'EdgeColor','None','FaceAlpha',1);
bar([3],nanmean(plot_data(:,3)),'FaceColor',[0.25,0.25,0.25],'EdgeColor','None','FaceAlpha',1);
errorbar([1],nanmean(plot_data(:,1)),(nanstd(plot_data(:,1))./sqrt(sum(~isnan(plot_data(:,1))))),'LineWidth',1,'Color',[0.75,0.75,0.75,1.0],'CapSize',0,'LineStyle','none')
errorbar([2],nanmean(plot_data(:,2)),(nanstd(plot_data(:,2))./sqrt(sum(~isnan(plot_data(:,2))))),'LineWidth',1,'Color',[0.5,0.5,0.5,1.0],'CapSize',0,'LineStyle','none')
errorbar([3],nanmean(plot_data(:,3)),(nanstd(plot_data(:,3))./sqrt(sum(~isnan(plot_data(:,3))))),'LineWidth',1,'Color',[0.25,0.25,0.25,1.0],'CapSize',0,'LineStyle','none')
if p_value(1) < 0.001
    title(['*** Left ', stim_choice_delay , ' weights'])
elseif p_value(1) < 0.01
    title(['** Left ', stim_choice_delay , ' weights'])
elseif p_value(1) < 0.05
    title(['* Left ', stim_choice_delay , ' weights'])
else
    title(['Left ', stim_choice_delay , ' weights'])
end
box(gca,'off')
ylim([0 0.05])
yticks([0 0.05])
ylabel('Weights')
xticks([])

% Plot right
subplot(1,2,2)
plot_data = animal_right_coeff;
bar([1],nanmean(plot_data(:,1)),'FaceColor',[0.75,0.75,0.75],'EdgeColor','None','FaceAlpha',1); hold on
bar([2],nanmean(plot_data(:,2)),'FaceColor',[0.5,0.5,0.5],'EdgeColor','None','FaceAlpha',1);
bar([3],nanmean(plot_data(:,3)),'FaceColor',[0.25,0.25,0.25],'EdgeColor','None','FaceAlpha',1);
errorbar([1],nanmean(plot_data(:,1)),(nanstd(plot_data(:,1))./sqrt(sum(~isnan(plot_data(:,1))))),'LineWidth',1,'Color',[0.75,0.75,0.75,1.0],'CapSize',0,'LineStyle','none')
errorbar([2],nanmean(plot_data(:,2)),(nanstd(plot_data(:,2))./sqrt(sum(~isnan(plot_data(:,2))))),'LineWidth',1,'Color',[0.5,0.5,0.5,1.0],'CapSize',0,'LineStyle','none')
errorbar([3],nanmean(plot_data(:,3)),(nanstd(plot_data(:,3))./sqrt(sum(~isnan(plot_data(:,3))))),'LineWidth',1,'Color',[0.25,0.25,0.25,1.0],'CapSize',0,'LineStyle','none')
if p_value(2) < 0.001
    title(['*** Right ', stim_choice_delay , ' weights'])
elseif p_value(2) < 0.01
    title(['** Right ', stim_choice_delay , ' weights'])
elseif p_value(2) < 0.05
    title(['* Right ', stim_choice_delay , ' weights'])
else
    title(['Right ', stim_choice_delay , ' weights'])
end
box(gca,'off')
ylim([0 0.05])
yticks([0 0.05])
ylabel('Weights')
xticks([])
end