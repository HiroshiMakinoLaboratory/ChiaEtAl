%% Get enrichment index across all regions, and across sessions.
function enrichment = get_enrichment_idx(parameters, data, coupling_sig_cells)
task_variable_cell_idx{1} = data.cell_idx.sample_left_union_sig_cells;
task_variable_cell_idx{2} = data.cell_idx.delay_left_union_sig_cells;
task_variable_cell_idx{3} = data.cell_idx.sample_right_union_sig_cells;
task_variable_cell_idx{4} = data.cell_idx.delay_right_union_sig_cells;
all_p_value = data.all_p_value;
total_num_animals = parameters.total_num_animals;
all_region_labels = data.all_region_labels;
total_num_session = 3; % Naive, inter, expert

for hide = 1:1 % Run and plot
    clear same_task_proportion stcr regional_same_task_proportion  mean_regional_stcr mean_stcr stcr len_task_related_coup_sig_cell stcr concat_stcr concat_num_task_related_coup_sig_cell concat_num_all_coup_sig_cell
    for k = 1:4
        for animal_num = 1: total_num_animals
            for session_num = 1: total_num_session
                same_task_proportion{animal_num}(:,session_num) = nan(length(all_p_value{animal_num}{session_num}) , 1);
                
                enrichment_index{k}{animal_num}{session_num} = nan(length(all_p_value{animal_num}{session_num}) , 1);
                len_task_related_coup_sig_cell{k}{animal_num}{session_num} = nan(length(all_p_value{animal_num}{session_num}) , 1);
                num_task_related_coup_sig_cell{k}{animal_num}{session_num} = nan(length(all_p_value{animal_num}{session_num}) , 1);
                num_all_coup_sig_cell{k}{animal_num}{session_num} = nan(length(all_p_value{animal_num}{session_num}) , 1);
                
                for cell_num = task_variable_cell_idx{k}{animal_num}{session_num}
                    
                    % Calculate enrichment index
                    if length(task_variable_cell_idx{k}{animal_num}{session_num}) >0
                        task_related_coup_sig_cell = (intersect(coupling_sig_cells{animal_num}{session_num}{cell_num}, task_variable_cell_idx{k}{animal_num}{session_num}));
                        
                        % If there are sig coupled cells, but no sig coupled cells to the specific task then STCR = 0
                        if (~isempty(task_related_coup_sig_cell) & sum(coupling_sig_cells{animal_num}{session_num}{cell_num}) ~=0)
                            
                            non_task_related_coup_sig_cell = setdiff(coupling_sig_cells{animal_num}{session_num}{cell_num},task_related_coup_sig_cell);
                            
                            len_task_related_coup_sig_cell{k}{animal_num}{session_num}(cell_num) = length((intersect(coupling_sig_cells{animal_num}{session_num}{cell_num}, task_variable_cell_idx{k}{animal_num}{session_num}))) / length(task_variable_cell_idx{k}{animal_num}{session_num});
                            num_task_related_coup_sig_cell{k}{animal_num}{session_num}(cell_num) = length(task_related_coup_sig_cell);
                            num_all_coup_sig_cell{k}{animal_num}{session_num}(cell_num) = length(coupling_sig_cells{animal_num}{session_num}{cell_num});
                            
                            enrichment_index{k}{animal_num}{session_num}(cell_num) = num_task_related_coup_sig_cell{k}{animal_num}{session_num}(cell_num) / num_all_coup_sig_cell{k}{animal_num}{session_num}(cell_num)';
                            
                            % If there are NO sig coupled cells then STCR = 0;
                        elseif sum(coupling_sig_cells{animal_num}{session_num}{cell_num}) == 0
                            enrichment_index{k}{animal_num}{session_num}(cell_num) = 0;
                        else
                            enrichment_index{k}{animal_num}{session_num}(cell_num) = 0;
                        end
                        clear task_related_coup_sig_cell non_task_related_coup_sig_cell
                    end
                end
            end
        end
        
        % Concat sessions across animals
        clear concat_all_region_labels
        for session_num = 1:total_num_session
            enrich_idx{k}{session_num} = [];
            concat_len_task_related_sig_cell{k}{session_num} = [];
            concat_all_region_labels{session_num} = [];
            concat_num_task_related_coup_sig_cell{k}{session_num} = [];
            concat_num_all_coup_sig_cell{k}{session_num} = [];
            for animal_num = 1:total_num_animals
                enrich_idx{k}{session_num} = [enrich_idx{k}{session_num} enrichment_index{k}{animal_num}{session_num}'];
                concat_len_task_related_sig_cell{k}{session_num} = [concat_len_task_related_sig_cell{k}{session_num} len_task_related_coup_sig_cell{k}{animal_num}{session_num}'];
                
                concat_all_region_labels{session_num} = [concat_all_region_labels{session_num} all_region_labels{animal_num}{session_num}'];
            end
        end
    end
end

% Averaging across animals
for animal_num = 1:13
    for session_num = 1:3
        sensory_concat_animal_stcr{animal_num}(:,session_num) =  [enrichment_index{1}{animal_num}{session_num}' enrichment_index{3}{animal_num}{session_num}'];
        choice_concat_animal_stcr{animal_num}(:,session_num) =  [enrichment_index{2}{animal_num}{session_num}' enrichment_index{4}{animal_num}{session_num}'];
    end
    mean_sensory_concat_animal_stcr(animal_num,:) = nanmean(sensory_concat_animal_stcr{animal_num},1);
    mean_choice_concat_animal_stcr(animal_num,:) = nanmean(choice_concat_animal_stcr{animal_num},1);
end

enrichment.averaged_choice = mean_choice_concat_animal_stcr;
enrichment.averaged_sensory = mean_sensory_concat_animal_stcr;
enrichment.enrich_idx = enrich_idx;

task_variable_cell_idx{1} = data.cell_idx.sample_left_union_sig_cells;
task_variable_cell_idx{2} = data.cell_idx.delay_left_union_sig_cells;
task_variable_cell_idx{3} = data.cell_idx.sample_right_union_sig_cells;
task_variable_cell_idx{4} = data.cell_idx.delay_right_union_sig_cells;

enrichment.sample.left_sample_enrich_idx = enrichment_index{1};
enrichment.sample.right_sample_enrich_idx = enrichment_index{3};
enrichment.choice.left_choice_enrich_idx = enrichment_index{2};
enrichment.choice.right_choice_enrich_idx = enrichment_index{4};
end
