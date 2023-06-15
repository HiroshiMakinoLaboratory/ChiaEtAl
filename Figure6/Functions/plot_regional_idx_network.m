%% Plot convergence index network
function plot_regional_idx_network(parameters,convergence_idx)
data1 = convergence_idx; % Initialize dataset

%% Intialize plotting hyperparameters
load('cortical_area_boundaries.mat')
bregma = [540,570];
bregma = bregma.*10;
region = 8;
region_coordinate_temp = [290,435; 350,270; 450,370; 465,500; 490,280; 630,200; 640,490; 690,340]; % Actual ALM: 380 shifted to 435 
region_coordinate_temp  = region_coordinate_temp*10;
region_coordinate_middle(:,1) = region_coordinate_temp(:,2) - bregma(2); % Swap X with Y
region_coordinate_middle(:,2) = region_coordinate_temp(:,1) - bregma(1);
region_coordinate_temp = region_coordinate_temp - bregma;
region_coordinate = region_coordinate_temp - bregma;
region_coordinate = (region_coordinate + 250)./10;

figure('Position',[200,100,200,200],'Color','w')
hold on;
% Plot gray background of cortex
for area_num = 1:size(cortical_area_boundaries)
    for hemi_num = 1:size(cortical_area_boundaries{area_num})
        plot(cortical_area_boundaries{area_num}{hemi_num}(:,1) - (bregma(1)./10),cortical_area_boundaries{area_num}{hemi_num}(:,2) - (bregma(2)./10),'LineWidth',1,'Color',[1,1,1,0.2])
        set(gca,'Color','k')
    end
end

% Plotting parameters
for hide = 1:1
    set(gca, 'YDir','reverse')
    legend off
    % plot bregma point
    hold on;
    ax = gca;
    set(ax,'FontSize',16)
    set(gca,'XDir','reverse');
    camroll(-270)
    xlim([-275 225]);
    ylim([-425 75]);
    box(ax,'on')
    set(ax,'XTick',[])
    set(ax,'YTick',[])
    hold on;
end

%% Plot network
% Make coolwarm colormap for lines
cMap = redblue(882*2);
cMap = cMap(882:1764,:); % Get red only colormap
colormap(cMap);
circle_multiplication_num = 14; % Arbitrary number used to scale the size of the circle
min_value= 0; max_value = 0.15;
multiplication_num = (882/max_value)-1; % Arbitrary number used to scale the colour and size of the arrows
region_idx = [1:region];

% Perform bootstrap to get significant increase in enrichment idex
num_shuffles = 1000; clear diff_p_value mean_diff concat_mean_diff
for region_num = 1:region
    for other_region_num = 1:region
        clear plot_data_expert plot_data_naive nan_idx remove_idx valid_idx
        % Cell x 1 vector (each value is enrich. index)
        plot_data_expert = data1{3}{region_num,other_region_num};
        plot_data_naive = data1{1}{region_num,other_region_num};
        plot_data_interm = data1{2}{region_num,other_region_num};
        
        % Get nan cells
        nan_idx = intersect(find(~isnan(plot_data_expert)), find(~isnan(plot_data_naive)));
        plot_data_expert = plot_data_expert(nan_idx);
        plot_data_naive = plot_data_naive(nan_idx);
                
        % Remove neurons that are consistently 0 across sessions
        remove_idx = intersect(intersect(find(plot_data_expert == 0),find(plot_data_naive == 0)),find(plot_data_interm==0) );
        valid_idx = [1:size(nan_idx,2)];
        valid_idx(remove_idx) = [];
        plot_data_expert = plot_data_expert(valid_idx);
        plot_data_naive = plot_data_naive(valid_idx);
        
        % Get averaged diff
        mean_diff(region_num,other_region_num) = mean(plot_data_expert -  plot_data_naive);
        
        rng(2022); clear shuffled_diff_mean
        for shuffle_num = 1:num_shuffles
            rand_idx = randsample(1:length(plot_data_expert),length(plot_data_expert),'true');
            shuffled_diff_mean(shuffle_num) = mean(plot_data_expert(rand_idx)-plot_data_naive(rand_idx));
        end
        diff_p_value(region_num,other_region_num) = sum(shuffled_diff_mean<0)/num_shuffles;
    end
end

% Get binary mask using significant p value
% Locate regions with sig. p values
p_value = diff_p_value;
p_value_a = (p_value < 0.05 & p_value >= 0.01) * 1/3;
p_value_b = (p_value < 0.01 & p_value >= 0.001) * 2/3;
p_value_c = (p_value < 0.001) * 3/3;
p_value_d = ones(length(p_value),length(p_value));
p_value_d(p_value >= 0.05) = nan;
p_value_combined_matrix = p_value_a+p_value_b+p_value_c .*p_value_d;
binary_mask = ~isnan(p_value_combined_matrix);

% Plot arrows with scaling factor
for other_region_num = region_idx
    for region_num = region_idx
        if other_region_num ~= region_num && binary_mask(region_num,other_region_num)
            scaling_factor = mean_diff(region_num,other_region_num);
            
            % Halfway coord
            source_coord = region_coordinate_temp(other_region_num,:)/10+25;
            target_coord = region_coordinate_temp(region_num,:)/10+25;
            
            % Calculate the correct point
            % Get the angle
            angle_titer = atand((source_coord(2)-target_coord(2)) / (source_coord(1)-target_coord(1)));
            h = sqrt( (source_coord(1)-target_coord(1))^2 + (source_coord(2)-target_coord(2))^2) - (circle_multiplication_num);
            o = sind((angle_titer))*h;
            a = cosd((angle_titer))*h;
            
            if other_region_num > region_num
                a = -a;
                o = -o;
            end
            
            adjusted_target_coord = [(source_coord(1) + a) (source_coord(2) + o)];
            
            if region_num ~= other_region_num
                if scaling_factor > max_value
                    scaling_factor = max_value;
                elseif scaling_factor < min_value
                    scaling_factor = min_value;
                end
                arrow( [source_coord(1), source_coord(2)], [adjusted_target_coord(1), adjusted_target_coord(2)],'Color',cMap( (ceil((scaling_factor*multiplication_num))) ,:),'length',6,'Width',round(abs(scaling_factor*5),2),'BaseAngle',90,'TipAngle',30)
            end
        end
    end
end

hold on;

% Plot Circle with scaling factor
diagonal_mask = diag(binary_mask);
for i = region_idx
    scaling_factor = mean_diff(i,i);
    if scaling_factor > max_value
        scaling_factor = max_value;
    elseif scaling_factor < min_value
        scaling_factor = min_value;
    end
    
    if diagonal_mask(i) == 0
        plot(region_coordinate_temp(i,1)/10+25,region_coordinate_temp(i,2)/10+25, 'o', 'MarkerSize',7.5,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor','k','LineWidth',1)
    else 
        plot(region_coordinate_temp(i,1)/10+25,region_coordinate_temp(i,2)/10+25, 'o', 'MarkerSize',7.5,'MarkerEdgeColor',cMap( (ceil((scaling_factor*multiplication_num))) ,:),'MarkerFaceColor',cMap( (ceil((scaling_factor*multiplication_num))) ,:),'LineWidth',1)
    end
end
axis square
% Add colorbar
c = colorbar;
c.Ticks = [0 1];
c.TickLabels = [{num2str(min_value) num2str(max_value)}];
a =  c.Position;
set(c,'Position',[a(1)+0.22 a(2) 0.04 0.82])
end