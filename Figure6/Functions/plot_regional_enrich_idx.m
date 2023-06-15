%% Plot regional enrichment index in heatmap
function plot_regional_enrich_idx(parameters, data1)
new_region_idx = parameters.new_region_idx;
region = 8; % 8 regions

% Perform bootstrap to get significant increase in enrichment idex
rng(2021)
num_shuffles = 1000; clear diff_p_value mean_diff concat_mean_diff
for region_num = 1:region
    for other_region_num = 1:region
        clear plot_data_expert plot_data_naive nan_idx remove_idx valid_idx
        % Cell x 1 vector (each value is enrich. index)
        plot_data_expert = data1{3}{region_num,other_region_num};
        plot_data_naive = data1{1}{region_num,other_region_num};
        plot_data_interm = data1{2}{region_num,other_region_num};
        
        % Get nan cells
        nan_idx = intersect(find(~isnan(plot_data_expert)),find(~isnan(plot_data_naive)));
        plot_data_expert = plot_data_expert(nan_idx);
        plot_data_naive = plot_data_naive(nan_idx);
        
        % Remove neurons that are consistently 0 across sessions
        remove_idx = intersect(intersect(find(plot_data_expert == 0),find(plot_data_naive == 0)),find(plot_data_interm==0));
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

% Get binary mask based on p_value
for hide = 1:1
    p_value = diff_p_value;
    c_axis = [0 1];
    % Locate regions with sig. p values
    p_value_a = (p_value < 0.05 & p_value >= 0.01) * 1/3;
    p_value_b = (p_value < 0.01 & p_value >= 0.001) * 2/3;
    p_value_c = (p_value < 0.001) * 3/3;
    p_value_d = ones(length(p_value),length(p_value));
    p_value_d(p_value >= 0.05) = nan;
    p_value_combined_matrix = p_value_a+p_value_b+p_value_c .*p_value_d;
    binary_mask = ~isnan(p_value_combined_matrix);
end

% Plot p value mask * increase in enrichment index
figure('Position',[200,200,150,150],'Color','black','DefaultAxesFontSize',14);
c_axis = [0 0.15];
for hide = 1:1
    cMap = redblue(100);
    cMap = cMap(51:100,:);
    plot_data = mean_diff(new_region_idx,new_region_idx);
    binary_mask_new = binary_mask(new_region_idx,new_region_idx);
    plot_data = plot_data.*binary_mask_new;
    plot_data(binary_mask_new ==0) = nan;
    imagesc(plot_data,'AlphaData',~isnan(plot_data))
    colormap(cMap)
    caxis(c_axis)
    xticks([]);yticks([])
    axis equal
    axis off
    c = colorbar;
    c.Color = [1 1 1];
end

% Adam Auton (2023). Red Blue Colormap (https://www.mathworks.com/matlabcentral/fileexchange/25536-red-blue-colormap), MATLAB Central File Exchange. Retrieved February 27, 2023.
    function c = redblue(m)
        %REDBLUE    Shades of red and blue color map
        %   REDBLUE(M), is an M-by-3 matrix that defines a colormap.
        %   The colors begin with bright blue, range through shades of
        %   blue to white, and then through shades of red to bright red.
        %   REDBLUE, by itself, is the same length as the current figure's
        %   colormap. If no figure exists, MATLAB creates one.
        %
        %   For example, to reset the colormap of the current figure:
        %
        %             colormap(redblue)
        %
        %   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG,
        %   COLORMAP, RGBPLOT.
        
        %   Adam Auton, 9th October 2009
        
        if nargin < 1, m = size(get(gcf,'colormap'),1); end
        
        if (mod(m,2) == 0)
            % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
            m1 = m*0.5;
            r = (0:m1-1)'/max(m1-1,1);
            g = r;
            r = [r; ones(m1,1)];
            g = [g; flipud(g)];
            b = flipud(r);
        else
            % From [0 0 1] to [1 1 1] to [1 0 0];
            m1 = floor(m*0.5);
            r = (0:m1-1)'/max(m1,1);
            g = r;
            r = [r; ones(m1+1,1)];
            g = [g; 1; flipud(g)];
            b = flipud(r);
        end
        
        c = [r g b];
        
    end

end