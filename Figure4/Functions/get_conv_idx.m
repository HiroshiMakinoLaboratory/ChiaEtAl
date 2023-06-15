%% Plot regional convergence index
function [convergence_idx] = get_conv_idx(data, coupling_sig_cells)
% Initialise data and parameters
total_num_animals = 13;
all_region_labels = data.all_region_labels;
task_variable_sig_cell = data.cell_idx.all_valid_cells; % Use valid cells in this analysis. Cells with significantly predicted activity from GLM.
total_num_session = 3; % Naive, intermediate and expert
region = 8; % 8 cortical regions

for animal_num = 1:total_num_animals
    for session_num = 1:total_num_session
        for region_num = 1:region
            % Check if there are cells to analyse. If there are no task_specific cells in the region. Value will be nan
            if (length(intersect(find(all_region_labels{animal_num}{session_num} == region_num),task_variable_sig_cell{animal_num}{session_num})') ~= 0)
                % Loop over each cell
                clear cell_idx
                for i = 1:length(intersect(find(all_region_labels{animal_num}{session_num} == region_num),task_variable_sig_cell{animal_num}{session_num})')
                    cell_idx = intersect(find(all_region_labels{animal_num}{session_num} == region_num),task_variable_sig_cell{animal_num}{session_num});
                    cell_num = cell_idx(i);
                    
                    %% Get intra-regional conv. index
                    % Target cell index is the task variable type for predictor cells
                    if ~isempty(intersect(find(all_region_labels{animal_num}{session_num} == region_num),intersect(coupling_sig_cells{animal_num}{session_num}{cell_num},task_variable_sig_cell{animal_num}{session_num})))
                        % Find cells within the same region that arealso the same-task behavioural variable. I.e similar to STCR
                        intra_region_cell_idx{animal_num}{session_num}{region_num}{i} = intersect(find(all_region_labels{animal_num}{session_num} == region_num),intersect(coupling_sig_cells{animal_num}{session_num}{cell_num},task_variable_sig_cell{animal_num}{session_num}));
                        num_all_coup_sig_cell= length(intersect(find(all_region_labels{animal_num}{session_num} == region_num)',coupling_sig_cells{animal_num}{session_num}{cell_num}));
                        
                        % Regional Conv. Idex -> Number of sig coupled cells in a predictor region / total num of cells in the predictor region
                        intra_region_idx{animal_num}{session_num}{region_num}(i) = num_all_coup_sig_cell / length(find(all_region_labels{animal_num}{session_num} == region_num)) ;
                    else
                        intra_region_idx{animal_num}{session_num}{region_num}(i) = 0;
                    end
                    
                    %% Get inter-regional conv. index
                    % Calculate other region - coupling e.g. Region 1-2, 1-3, 1-4.. etc
                    other_region_idx = [1:region]; other_region_idx(region_num) = [];
                    
                    for other_region_num = other_region_idx
                        if ~isempty(intersect(find(all_region_labels{animal_num}{session_num} == other_region_num),intersect(coupling_sig_cells{animal_num}{session_num}{cell_num},task_variable_sig_cell{animal_num}{session_num})))
                            other_region_cell_idx{animal_num}{session_num}{region_num}{other_region_num}{i} = intersect(find(all_region_labels{animal_num}{session_num} == other_region_num),intersect(coupling_sig_cells{animal_num}{session_num}{cell_num},task_variable_sig_cell{animal_num}{session_num}));
                            
                            num_task_related_coup_sig_cell= length(other_region_cell_idx{animal_num}{session_num}{region_num}{other_region_num}{i} );
                            num_all_coup_sig_cell= length(intersect(find(all_region_labels{animal_num}{session_num} == other_region_num)',coupling_sig_cells{animal_num}{session_num}{cell_num}));
                            
                            % Regional Conv. Idex -> Number of sig coupled cells in a predictor region / total num of cells in the predictor region
                            other_region_cell_prop{animal_num}{session_num}{region_num,other_region_num}(i) = num_all_coup_sig_cell / length(find(all_region_labels{animal_num}{session_num} == other_region_num));
                            
                            clear num_all_coup_sig_cell num_task_related_coup_sig_cell
                        else
                            % If there are no intersect between region, task_variable and coupling, proportion will be 0
                            other_region_cell_prop{animal_num}{session_num}{region_num,other_region_num}(i) = 0;
                        end
                    end
                    clear other_region_idx non_task_variable_idx non_intra_region_cell_idx non_intra_region_idx non_other_region_cell_idx non_other_region_idx
                end
                
            else
                % If there are no task specific cells in the regions. All the values will be nan, including the same region and the other regions.
                other_region_idx = [1:region]; other_region_idx(region_num) = [];
                
                intra_region_cell_idx{animal_num}{session_num}{region_num} = nan;
                intra_region_idx{animal_num}{session_num}{region_num} = nan;
                for other_region_num = other_region_idx
                    other_region_cell_idx{animal_num}{session_num}{region_num}{other_region_num} = nan;
                    other_region_cell_prop{animal_num}{session_num}{region_num,other_region_num} = nan;
                end
            end
        end
    end
end

% Concat into its session and regional variables
for session_num = 1:total_num_session
    for region_num = 1:region
        concat_regional_prop_intra{session_num}{region_num} = [];
        
        for animal_num = 1:total_num_animals
            concat_regional_prop_intra{session_num}{region_num} = [concat_regional_prop_intra{session_num}{region_num} intra_region_idx{animal_num}{session_num}{region_num}];
        end
    end
end

% Concat other regions prop into its session and regional variables
for session_num = 1:total_num_session
    for region_num = 1:region
        other_region_idx = [1:region]; other_region_idx(region_num) = [];
        concat_regional_prop_other{session_num}{region_num,region_num} = concat_regional_prop_intra{session_num}{region_num};
        
        for other_region_num = other_region_idx
            concat_regional_prop_other{session_num}{region_num,other_region_num} = [];
            for animal_num = 1:total_num_animals
                concat_regional_prop_other{session_num}{region_num,other_region_num} = [concat_regional_prop_other{session_num}{region_num,other_region_num} other_region_cell_prop{animal_num}{session_num}{region_num,other_region_num}];
            end
            matrix_mean_concat_regional_prop_other{session_num}(region_num,other_region_num) = nanmean(concat_regional_prop_other{session_num}{region_num,other_region_num});
        end
        matrix_mean_concat_regional_prop_other{session_num}(region_num,region_num) = nanmean(concat_regional_prop_other{session_num}{region_num,region_num});
        
        clear other_region_idx
    end
end

% Concat across animals
for session_num = 1:3
    for animal_num = 1:13
        for region_num = 1:8
            for other_region_num = 1:8
                if region_num == other_region_num
                    concat_animals_conv_idx{session_num}(region_num,other_region_num,animal_num) = nanmean(intra_region_idx{animal_num}{session_num}{region_num});
                else
                    concat_animals_conv_idx{session_num}(region_num,other_region_num,animal_num) = nanmean(other_region_cell_prop{animal_num}{session_num}{region_num,other_region_num});
                end
            end
        end
    end
end

% Analyse and reshape intra vs inter region coupling
for session_num = 1:total_num_session
    for region_num = 1:region
        
        % Concat across regions, so each cell will have an other region_conv_idx
        other_region_vector = setdiff([1:region],region_num);
        for  i = 1:length(other_region_vector)
            other_region_num = other_region_vector(i);
            reshaped_vector(:,:,i) = concat_regional_prop_other{session_num}{region_num,other_region_num};
        end
        
        % get mean of inter regional conv idx
        concat_all_other_region_conv_idx{session_num}{region_num} = squeeze(nanmean(reshaped_vector,3));
        
        intra_only{session_num}{region_num} = concat_regional_prop_other{session_num}{region_num,region_num} ;
        inter_only{session_num}{region_num} = concat_all_other_region_conv_idx{session_num}{region_num};
        
        clear other_region_vector temp reshaped_vector
    end
    
    inter_intra.mean_intra_only{session_num} = cellfun(@nanmean,intra_only{session_num});
    inter_intra.sem_intra_only{session_num} = cellfun(@nanstd,intra_only{session_num}) ./ sqrt(cellfun(@length,intra_only{session_num}));
    
    inter_intra.mean_inter_only{session_num} = cellfun(@nanmean,inter_only{session_num});
    inter_intra.sem_inter_only{session_num} = cellfun(@nanstd,inter_only{session_num}) ./ sqrt(cellfun(@length,inter_only{session_num}));
end

% Store variables in struct
convergence_idx.raw.intra_only = intra_only;
convergence_idx.raw.inter_only = inter_only;
convergence_idx.matrix_mean_concat_regional_prop_other =  matrix_mean_concat_regional_prop_other;
convergence_idx.inter_intra.mean_intra_only =  inter_intra.mean_intra_only;
convergence_idx.inter_intra.sem_intra_only =  inter_intra.sem_intra_only;
convergence_idx.inter_intra.mean_inter_only =  inter_intra.mean_inter_only;
convergence_idx.inter_intra.sem_inter_only =  inter_intra.sem_inter_only;
convergence_idx.regional.concat_animals_conv_idx = concat_animals_conv_idx;

%% Get global convergence index (across all regions)

all_region_labels = data.all_region_labels;
all_valid_cells = data.cell_idx.all_valid_cells;
for animal_num = 1:size(all_region_labels,2)
    for session_num = 1:3
        conv_idxs{animal_num}{session_num} = nan(1,size(coupling_sig_cells{animal_num}{session_num},2));
        for cell_num = all_valid_cells{animal_num}{session_num}
            if ~isempty(coupling_sig_cells{animal_num}{session_num}{1})
                
                conv_idxs{animal_num}{session_num}(cell_num) = length(coupling_sig_cells{animal_num}{session_num}{cell_num}) ./ size(coupling_sig_cells{animal_num}{session_num},2);
            else
                conv_idxs{animal_num}{session_num}(cell_num) = nan;
            end
        end
        % Get average across animals
        mean_conv_idxs(animal_num,session_num) =  nanmean(conv_idxs{animal_num}{session_num});
    end
end

% Concat across all animals
for session_num = 1:3
    concat_conv_idx{session_num}  = [];
    for animal_num = 1:size(all_region_labels,2)
        concat_conv_idx{session_num} = [concat_conv_idx{session_num} conv_idxs{animal_num}{session_num}];
    end
end

convergence_idx.global.concat_conv_idx = concat_conv_idx;
convergence_idx.global.conv_idxs = conv_idxs;
convergence_idx.global.mean_conv_idxs = mean_conv_idxs;
end