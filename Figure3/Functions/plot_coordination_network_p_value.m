%% Plot coordination network after testing for significance
function p_value_matrix = plot_coordination_network_p_value(parameters,coordination_data, sensory_or_choice, both_or_same, session_1, session_2)

%% Initialise data based on parameters
if strcmp(both_or_same,'same')
    data1 = coordination_data.(['trial_by_trial_corr_',sensory_or_choice,'_left_',session_1,'_concat']);
    data2 = coordination_data.(['trial_by_trial_corr_',sensory_or_choice,'_right_',session_1,'_concat']);
    data3 = coordination_data.(['trial_by_trial_corr_',sensory_or_choice,'_left_',session_2,'_concat']);
    data4 = coordination_data.(['trial_by_trial_corr_',sensory_or_choice,'_right_',session_2,'_concat']);
    
    % Concat left and right data
    data_naive =  cat(1,data1,data2);
    data_expert = cat(1,data3,data4);
else
    data_naive = coordination_data.(['trial_by_trial_corr_',sensory_or_choice,'_',both_or_same,'_',session_1,'_concat']);
    data_expert = coordination_data.(['trial_by_trial_corr_',sensory_or_choice,'_',both_or_same,'_',session_2,'_concat']);
end

% Do bootstrap to get p values
for region_num1 = 1:8
    for region_num2 = 1:8
        clear naive_temp naive expert_temp expert naive_shuffle expert_shuffle
        
        naive_temp = data_naive(:,region_num1,region_num2);
        naive = naive_temp(~isnan(naive_temp));
        expert_temp = data_expert(:,region_num1,region_num2);
        expert = expert_temp(~isnan(expert_temp));
        
        for shuffle_num = 1:1000
            for session_num = 1:length(naive)
                naive_shuffle(shuffle_num,session_num) = naive(randi(length(naive)));
            end
            for session_num = 1:length(expert)
                expert_shuffle(shuffle_num,session_num) = expert(randi(length(expert)));
            end
        end
        p_value_matrix(region_num1,region_num2) = sum(mean(naive_shuffle,2) > mean(expert_shuffle,2))./1000;
    end
end

% Remove intra-regional coordination
for region_num = 1:8
    for region_num = 1:8
        p_value_matrix(region_num,region_num) = 1; % Change to 1
    end
end

% Bin the p-values to 0.05, 0.01 and 0.001
p_value_bin = double(p_value_matrix < 0.001) + double(p_value_matrix < 0.01) + double(p_value_matrix < 0.05);
upper_p_value_bin = triu(p_value_bin,1);
p_value_bin = upper_p_value_bin + upper_p_value_bin';

for region_num = 1:8
    for region_num = 1:8
        p_value_bin(region_num,region_num) = 0;
    end
end

%% Plot in real coordinate.
% Intialise plotting hyperparameters
load('cortical_area_boundaries.mat') % Load cortical area map
bregma = [540,570]; % Fixed number based on cortical area map
bregma = bregma.*10;

% Get coordinates of all regions.
region_coordinate_temp = [290,435; 350,270; 450,370; 465,500; 490,280; 630,200; 640,490; 690,340]; % Actual ALM: 380 shifted to 435
region_coordinate_temp  = region_coordinate_temp.*10;
region_coordinate_middle(:,1) = region_coordinate_temp(:,2) - bregma(2); % Swap X with Y.
region_coordinate_middle(:,2) = region_coordinate_temp(:,1) - bregma(1);
region_coordinate = region_coordinate_temp - bregma;

% Plot.
figure('Position',[200,200,200,200],'Color','w')
hold on;
for area_num = 1:size(cortical_area_boundaries)
    for hemi_num = 1:size(cortical_area_boundaries{area_num})
        plot(cortical_area_boundaries{area_num}{hemi_num}(:,1) - (bregma(1)./10),cortical_area_boundaries{area_num}{hemi_num}(:,2) - (bregma(2)./10),'LineWidth',1,'Color',[1,1,1,0.2])
    end
end

% Plot parameters
xlim([-275 225]);
ylim([-425 75]);
box on
axis square
ax = gca;
ax.Color = 'k';
set(gca,'XDir','reverse');
set(gca,'YDir','reverse')
set(ax,'XTick',[])
set(ax,'YTick',[])
camroll(-270)
cMap = redblue(11);
cMap = cMap([7 9 11],:);
line_idx = p_value_bin;

region_coordinate = (region_coordinate + 250)./10;
for region_num1 = 1:8
    for region_num2 = 1:8
        
        new_line_idx = abs(line_idx(region_num1,region_num2));
        if p_value_bin(region_num1,region_num2) > 0 % Only display significant pairs
            line([region_coordinate(region_num1,1),region_coordinate(region_num2,1)],[region_coordinate(region_num1,2),region_coordinate(region_num2,2)],'LineWidth',new_line_idx./1.5,'Color',cMap(new_line_idx,:))
        end
    end
end

for region_num = 1:8
    plot(region_coordinate(region_num,1),region_coordinate(region_num,2),'o','MarkerSize',64./8,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor','None')
end
end