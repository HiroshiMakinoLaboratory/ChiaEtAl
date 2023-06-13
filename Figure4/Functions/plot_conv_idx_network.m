%% Plot convergence index network
function plot_conv_idx_network(parameters,convergence_idx)
%% Intialise plotting hyperparameters
cd(parameters.root_folder)
load('cortical_area_boundaries.mat')

% Load cortical area boundaries mat file required to plot outline of cortical areas
bregma = [540,570];
bregma = bregma.*10;
region = 8;

% Coordinates for 8 regions
region_coordinate_temp = [290,435; 350,270; 450,370; 465,500; 490,280; 630,200; 640,490; 690,340]; % Actual ALM: 380 shifted to 435
region_coordinate_temp  = region_coordinate_temp*10;
region_coordinate_middle(:,1) = region_coordinate_temp(:,2) - bregma(2); % Swap X with Y
region_coordinate_middle(:,2) = region_coordinate_temp(:,1) - bregma(1);
region_coordinate_temp = region_coordinate_temp - bregma;

for session_num = 1:3
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
    cMap = colormap(magma(882)); % Arbitrary value is used to scale the colours
    circle_multiplication_num = 7; % Arbitrary number used to scale the size of the circle
    min_value= 0.05; max_value = 0.25;
    multiplication_num = (882/max_value)-1; % Arbitrary number used to scale the colour and size of the arrows
    region_idx = [1:region];
    % Input data to use
    data_to_use = convergence_idx.matrix_mean_concat_regional_prop_other{session_num};
    
    % Plot arrows with scaling factor
    for other_region_num = region_idx
        for region_num = region_idx
            
            if other_region_num ~= region_num
                scaling_factor = data_to_use(region_num,other_region_num);
                
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
                    
                    arrow( [source_coord(1), source_coord(2)], [adjusted_target_coord(1), adjusted_target_coord(2)],'Color',cMap( (ceil((scaling_factor*multiplication_num))) ,:),'length',6,'Width',abs(scaling_factor*11),'BaseAngle',90,'TipAngle',30)
                end
            end
        end
    end
    
    hold on;
    
    % Plot Circle with scaling factor
    scaling_vector = diag(data_to_use); % Chose only intra-regional numbers
    for i = region_idx
        scaling_num = scaling_vector(i);
        if scaling_num < min_value
            scaling_num = min_value
        elseif scaling_num >= max_value
            scaling_num = max_value;
        end
        plot(region_coordinate_temp(i,1)/10+25,region_coordinate_temp(i,2)/10+25, '.', 'MarkerSize',25,'Color',cMap(ceil(scaling_num*multiplication_num),:),'MarkerFaceColor',cMap(ceil(scaling_num*multiplication_num),:)  )
    end
    axis square
    
    if session_num == 1
        title('Naive')
    elseif session_num == 2
        title('Intermediate')
    elseif session_num == 3
        title('Expert')
    end
    
    % Add colorbar
    c = colorbar;
    c.Ticks = [0 1];
    c.TickLabels = [{num2str(min_value) num2str(max_value)}];
    a =  c.Position;
    set(c,'Position',[a(1)+0.22 a(2) 0.04 0.82])
end
end