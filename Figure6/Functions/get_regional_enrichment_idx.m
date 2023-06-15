%% Get regional enrichment idex within task-specific neurons
% i.e Delay cells in region 1 to delay cells in region 8
function regional_enrichment = get_regional_enrichment_idx(parameters,data, coupling_sig_cells)
task_variable_cell_idx{1} = data.cell_idx.sample_right_union_sig_cells;
task_variable_cell_idx{2} = data.cell_idx.delay_right_union_sig_cells;
region = 8; % Total number of regions used
total_num_animals = parameters.total_num_animals;
total_num_session = 3; % Naive, inter ,expert
all_region_labels = data.all_region_labels;

for hide = 1:1
    clear concat_regional_prop_other intra_region_idx other_region_cell_idx concat_regional_prop_intra other_region_cell_prop intra_region_matrix intra_region_cell_idx all_inter_region_cell_idx regional_inter_region_cell_idx inter_region_cell_idx intra_region_metric intra_region_cell_idx intra_region_matrix other_region_metric   intra_region_matrix
    for target_cell_idx = 1:2 % Target index e.g. delay left cells
        for k = 1:2 % Source index e.g. delay left cells
            intra_region_metric = [];  other_region_metric=[];
            for animal_num = 1:size(task_variable_cell_idx{1},2)
                
                for session_num = 1:3
                    for region_num = 1:region
                        
                        % Check if there are cells to analyse. If there are no task_specific cells in the region. Value will be nan
                        if (length(intersect(find(all_region_labels{animal_num}{session_num} == region_num),task_variable_cell_idx{k}{animal_num}{session_num})') ~= 0)
                            % Loop over each cell
                            clear cell_idx
                            % First initialise the cell
                            intra_region_idx{k}{animal_num}{session_num}{region_num} = nan(length(intersect(find(all_region_labels{animal_num}{session_num} == region_num),task_variable_cell_idx{k}{animal_num}{session_num})'),1)';
                            % Calculate other region - coupling e.g. Region 1-2, 1-3, 1-4.. etc
                            other_region_idx = [1:region]; other_region_idx(region_num) = [];
                            for other_region_num = other_region_idx
                                other_region_cell_prop{k}{animal_num}{session_num}{region_num,other_region_num} = nan(length(intersect(find(all_region_labels{animal_num}{session_num} == region_num),task_variable_cell_idx{k}{animal_num}{session_num})'),1)';
                            end
                            
                            for i = 1:length(intersect(find(all_region_labels{animal_num}{session_num} == region_num),task_variable_cell_idx{k}{animal_num}{session_num})')
                                cell_idx = intersect(find(all_region_labels{animal_num}{session_num} == region_num),task_variable_cell_idx{k}{animal_num}{session_num});
                                cell_num = cell_idx(i);
                                % Target cell index is the task variable type for predictor cells
                                % check if there are valid coupling cells first
                                if ~isnan(sum(coupling_sig_cells{animal_num}{session_num}{cell_num}))   | isempty(coupling_sig_cells{animal_num}{session_num}{cell_num})
                                    if ~isempty(intersect(find(all_region_labels{animal_num}{session_num} == region_num),intersect(coupling_sig_cells{animal_num}{session_num}{cell_num},task_variable_cell_idx{target_cell_idx}{animal_num}{session_num})))
                                        % Find cells within the same region that arealso the same-task behavioural variable. I.e similar to STCR
                                        intra_region_cell_idx{k}{animal_num}{session_num}{region_num}{i} = intersect(find(all_region_labels{animal_num}{session_num} == region_num),intersect(coupling_sig_cells{animal_num}{session_num}{cell_num},task_variable_cell_idx{target_cell_idx}{animal_num}{session_num}));
                                        inter_region_cell_idx{k}{animal_num}{session_num}{region_num}{i} = intersect(find(all_region_labels{animal_num}{session_num} ~= region_num),intersect(coupling_sig_cells{animal_num}{session_num}{cell_num},task_variable_cell_idx{target_cell_idx}{animal_num}{session_num}));
                                        
                                        num_task_related_coup_sig_cell = length(intra_region_cell_idx{k}{animal_num}{session_num}{region_num}{i});
                                        num_all_coup_sig_cell= length(intersect(find(all_region_labels{animal_num}{session_num} == region_num)',coupling_sig_cells{animal_num}{session_num}{cell_num}));
                                        
                                        % Get enrichment index by normalising with all sig coupled cells
                                        intra_region_idx{k}{animal_num}{session_num}{region_num}(i) = num_task_related_coup_sig_cell /num_all_coup_sig_cell;
                                        
                                        % If there are coupling predictors, but none of them interesect with the cell's predictors
                                    else
                                        intra_region_idx{k}{animal_num}{session_num}{region_num}(i) = 0;
                                    end
                                else
                                    intra_region_idx{k}{animal_num}{session_num}{region_num}(i) = nan;
                                end
                                
                                % Calculate other region - coupling e.g. Region 1-2, 1-3, 1-4.. etc
                                other_region_idx = [1:region]; other_region_idx(region_num) = [];
                                
                                for other_region_num = other_region_idx
                                    if ~isnan(sum(coupling_sig_cells{animal_num}{session_num}{cell_num})) | isempty(coupling_sig_cells{animal_num}{session_num}{cell_num})
                                        if ~isempty(intersect(find(all_region_labels{animal_num}{session_num} == other_region_num),intersect(coupling_sig_cells{animal_num}{session_num}{cell_num},task_variable_cell_idx{target_cell_idx}{animal_num}{session_num})))
                                            other_region_cell_idx{k}{animal_num}{session_num}{region_num}{other_region_num}{i} = intersect(find(all_region_labels{animal_num}{session_num} == other_region_num),intersect(coupling_sig_cells{animal_num}{session_num}{cell_num},task_variable_cell_idx{target_cell_idx}{animal_num}{session_num}));
                                            num_task_related_coup_sig_cell= length(other_region_cell_idx{k}{animal_num}{session_num}{region_num}{other_region_num}{i} );
                                            num_all_coup_sig_cell= length(intersect(find(all_region_labels{animal_num}{session_num} == other_region_num)',coupling_sig_cells{animal_num}{session_num}{cell_num}));
                                            other_region_cell_prop{k}{animal_num}{session_num}{region_num,other_region_num}(i) = num_task_related_coup_sig_cell / num_all_coup_sig_cell;
                                            
                                            clear a b c union_a_b_c
                                            clear num_all_coup_sig_cell num_task_related_coup_sig_cell
                                        else
                                            % If there are no intersect between region, task_variable and coupling, proportion will be 0
                                            other_region_cell_idx{k}{animal_num}{session_num}{region_num}{other_region_num}{i} = 0;
                                            other_region_cell_prop{k}{animal_num}{session_num}{region_num,other_region_num}(i) = 0;
                                        end
                                    else
                                        intra_region_idx{k}{animal_num}{session_num}{region_num}(i) = nan;
                                    end
                                end
                                clear other_region_idx non_task_variable_idx non_intra_region_cell_idx non_intra_region_idx non_other_region_cell_idx non_other_region_idx
                            end
                            
                        else
                            % If there are no task specific cells in the regions. All the values will be nan, including the same region and the other regions.
                            other_region_idx = [1:region]; other_region_idx(region_num) = [];
                            intra_region_cell_idx{k}{animal_num}{session_num}{region_num} = nan;
                            intra_region_idx{k}{animal_num}{session_num}{region_num} = nan;
                            for other_region_num = other_region_idx
                                other_region_cell_idx{k}{animal_num}{session_num}{region_num}{other_region_num} = nan;
                                other_region_cell_prop{k}{animal_num}{session_num}{region_num,other_region_num} = nan;
                            end
                        end
                    end
                end
            end
            
            % Concat into its session and regional variables
            for session_num = 1:total_num_session
                for region_num = 1:region
                    
                    % Initialise
                    concat_regional_prop_intra{k}{session_num}{region_num} = [];
                    for other_region_num = 1:region
                        for animal_num = 1:total_num_animals
                            concat_animal_regional_prop_other{target_cell_idx}{k}{session_num}{region_num,other_region_num}{animal_num} = [];
                        end
                    end
                    
                    animal_idx_enrich_idx{target_cell_idx}{k}{session_num}{region_num,region_num} = [];
                    
                    for animal_num = 1:total_num_animals
                        concat_regional_prop_intra{k}{session_num}{region_num} = [concat_regional_prop_intra{k}{session_num}{region_num} intra_region_idx{k}{animal_num}{session_num}{region_num}];
                        concat_animal_regional_prop_other{target_cell_idx}{k}{session_num}{region_num,region_num}{animal_num} = [concat_animal_regional_prop_other{target_cell_idx}{k}{session_num}{region_num,region_num}{animal_num}  intra_region_idx{k}{animal_num}{session_num}{region_num}];
                        animal_idx_enrich_idx{target_cell_idx}{k}{session_num}{region_num,region_num} =  [animal_idx_enrich_idx{target_cell_idx}{k}{session_num}{region_num,region_num}  repelem(animal_num,length(intra_region_idx{k}{animal_num}{session_num}{region_num}) )];
                    end
                end
            end
            
            % Concat other regions prop into its session and regional variables
            for session_num = 1:total_num_session
                for region_num = 1:region
                    other_region_idx = [1:region]; other_region_idx(region_num) = [];
                    concat_regional_prop_other{target_cell_idx}{k}{session_num}{region_num,region_num} = concat_regional_prop_intra{k}{session_num}{region_num};
                    
                    for other_region_num = other_region_idx
                        concat_regional_prop_other{target_cell_idx}{k}{session_num}{region_num,other_region_num} = [];
                        
                        % Manually initialise
                        animal_idx_enrich_idx{target_cell_idx}{k}{session_num}{region_num,other_region_num} = [];
                        for animal_num = 1:total_num_animals
                            concat_regional_prop_other{target_cell_idx}{k}{session_num}{region_num,other_region_num} = [concat_regional_prop_other{target_cell_idx}{k}{session_num}{region_num,other_region_num} other_region_cell_prop{k}{animal_num}{session_num}{region_num,other_region_num}];
                            animal_idx_enrich_idx{target_cell_idx}{k}{session_num}{region_num,other_region_num} =  [animal_idx_enrich_idx{target_cell_idx}{k}{session_num}{region_num,other_region_num}  repelem(animal_num,length(other_region_cell_prop{k}{animal_num}{session_num}{region_num,other_region_num}) )];
                        end
                    end
                    clear other_region_idx
                end
            end
        end
    end
end

regional_enrichment.delay_right = concat_regional_prop_other{2}{2};
regional_enrichment.stim_right = concat_regional_prop_other{1}{1};

regional_enrichment.animal_idx.delay_right = animal_idx_enrich_idx{2}{2};
regional_enrichment.animal_idx.stim_right = animal_idx_enrich_idx{1}{1};
end