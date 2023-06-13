%% Plot epoch-average choice activity
function plot_ablated_epoch_averaged_choice_mode(parameters, ablated_data)
session_num = 1;
fs_image = parameters.fs_image;
x_lim = [round(fs_image*2),round(fs_image*4)]; y_lim = [-2 2];
region_idx = [1 6 7 8]; % Regions ALM, vS1, RSC and PPC
without_resp_epoch = parameters.pre_resp_epoch;
stim_offset = parameters.sample_offset;

concat_original_right_choice_mode = ablated_data.choice_mode.concat_original_right_choice_mode;
concat_original_left_choice_mode = ablated_data.choice_mode.concat_original_left_choice_mode;
concat_control_right_choice_mode = ablated_data.choice_mode.concat_control_right_choice_mode;
concat_control_left_choice_mode = ablated_data.choice_mode.concat_control_left_choice_mode;

decode_choice_epoch = parameters.choice_epoch;

for region_num = region_idx
    for other_region_num = region_idx
        for cell_num_removal = [1]
            concat_diff_data_1=[];concat_diff_data_2=[];concat_diff_data_3=[];concat_diff_data_4=[];
            valid_animal_idx = find(abs(cellfun(@isempty,concat_original_right_choice_mode{session_num}{region_num,other_region_num}{cell_num_removal})-1));
            
            for animal_num = valid_animal_idx
                for shuffle_num = 1:size(concat_original_right_choice_mode{session_num}{region_num,other_region_num}{cell_num_removal}{animal_num},3)
                    original_data_right = concat_original_right_choice_mode{session_num}{region_num,other_region_num}{cell_num_removal}{animal_num}(:,without_resp_epoch,shuffle_num);
                    original_data_left = concat_original_left_choice_mode{session_num}{region_num,other_region_num}{cell_num_removal}{animal_num}(:,without_resp_epoch,shuffle_num);
                    
                    ablated_data_right = concat_control_right_choice_mode{session_num}{region_num,other_region_num}{cell_num_removal}{animal_num}(:,without_resp_epoch,shuffle_num);
                    ablated_data_left = concat_control_left_choice_mode{session_num}{region_num,other_region_num}{cell_num_removal}{animal_num}(:,without_resp_epoch,shuffle_num);
                    
                    % Average across trials
                    original_data_right = nanmean(original_data_right,1);
                    original_data_left = nanmean(original_data_left,1);
                    
                    ablated_data_right = nanmean(ablated_data_right,1);
                    ablated_data_left = nanmean(ablated_data_left,1);
                    
                    % Offset each data
                    offset_original_right  = nanmean(original_data_right(:,stim_offset),2);
                    offset_original_left  = nanmean(original_data_left(:,stim_offset),2);
                    
                    offset_ablated_right  = nanmean(ablated_data_right(:,stim_offset),2);
                    offset_ablated_left = nanmean(ablated_data_left(:,stim_offset),2);
                    
                    % adjust according to remove offset
                    original_data_right = original_data_right - offset_original_right;
                    original_data_left = original_data_left - offset_original_left;
                    
                    ablated_data_right = ablated_data_right - offset_ablated_right;
                    ablated_data_left = ablated_data_left - offset_ablated_left;
                    
                    % Max normalise
                    norm_original_data_right = original_data_right ./nanstd([original_data_right ablated_data_right   ]);
                    norm_original_data_left = original_data_left ./nanstd([original_data_right ablated_data_right   ]);
                    
                    norm_ablated_data_right = ablated_data_right ./nanstd([original_data_right ablated_data_right   ]);
                    norm_ablated_data_left = ablated_data_left ./nanstd([original_data_right ablated_data_right    ]);
                    
                    % Get the difference
                    diff_data_1(:,shuffle_num) = ((norm_original_data_right ))';
                    diff_data_2(:,shuffle_num) = ((norm_ablated_data_right ))';
                    diff_data_3(:,shuffle_num) = ((norm_original_data_left))';
                    diff_data_4(:,shuffle_num) = ((norm_ablated_data_left ))';
                    
                end
                mean_diff_data_1 = nanmean(diff_data_1,2);
                mean_diff_data_2 = nanmean(diff_data_2,2);
                mean_diff_data_3 = nanmean(diff_data_3,2);
                mean_diff_data_4 = nanmean(diff_data_4,2);
                
                concat_diff_data_1 = [concat_diff_data_1 mean_diff_data_1];
                concat_diff_data_2 = [concat_diff_data_2 mean_diff_data_2];
                concat_diff_data_3 = [concat_diff_data_3 mean_diff_data_3];
                concat_diff_data_4 = [concat_diff_data_4 mean_diff_data_4];
                
                clear diff_data_1 diff_data_2 diff_data_3
            end
            choice_concat_mean_1{region_num,other_region_num} =nanmean(concat_diff_data_1(decode_choice_epoch,:));
            choice_concat_mean{region_num,other_region_num} = nanmean(concat_diff_data_2(decode_choice_epoch,:)) - nanmean(concat_diff_data_1(decode_choice_epoch,:));
        end
    end
end

%% Plot bar plot for Target:ALM
figure('Position',[200,100,150,170],'Color','white','DefaultAxesFontSize',14);
counter = 1;
for other_region_num = [6 7 8]
    
    for region_num = 1
        hold on;
        a = choice_concat_mean{region_num,other_region_num} ;
        bar(counter,nanmean(a),'FaceColor','k','EdgeColor','None','FaceAlpha',1)
        errorbar(counter,nanmean(a),sqrt(nanstd(a)./length(a)),'LineWidth',1,'Color','k','CapSize',0,'LineStyle','none')
        counter = counter+1;
        
        n_num(other_region_num) = length(a);
    end
    xlim([0 4])
    ylim([-2 1])
    ylabel('Diff. Choice activity')
    xticks([1 2 3])
    xticklabels({'vS1','RSC','PPC'})
    set(gca,'XTickLabelRotation',90)
    title('Target: ALM')
end

%% Plot heatmap
clear plot_data

% Remove intra-regional ablations and remove nan regions
counter_b = 1;
for region_num = [1 6 7 8]
    counter_a = 1;
    for other_region_num = [1 6 7 8]
        plot_data(counter_b,counter_a) = nanmean(choice_concat_mean{region_num,other_region_num});
        counter_a = counter_a +1;
    end
    counter_b = counter_b+1;
end

% Do stats bootstrap to get significant increase in choice activity
for region_num = [1 6 7 8]
    for other_region_num = [1 6 7 8]
        for shuffle_num = 1:1000
            rand_idx =randsample(length(choice_concat_mean{region_num,other_region_num}) , length(choice_concat_mean{region_num,other_region_num}),'true');
            choice_shuffle(shuffle_num) = nanmean(choice_concat_mean{region_num,other_region_num}(rand_idx));
        end
        p_value(region_num,other_region_num) = sum(choice_shuffle>0)/1000;
    end
end
all_p_value = p_value([1 6 7 8],[1 6 7 8]);

% Find data with p value < 0.025, two sided test
neg_sig_data = plot_data .* (all_p_value < (0.05/2));
pos_sig_data = plot_data .* (all_p_value > ( (1-0.05/2)) );
sig_data = neg_sig_data+ pos_sig_data;
sig_data(find(sig_data ==0)) = nan;

% Plot heatmap of differences only sig ones > 0.025
cMap = redblue(100);
figure('Position',[200,200,150,150],'Color','black','DefaultAxesFontSize',14);
imagesc(sig_data,'AlphaData',~isnan(sig_data));
c_axis = [-1.0 1.0];
colormap(cMap)
caxis(c_axis)
xticks([]);yticks([]);
axis equal
axis off
cbar = colorbar;
cbar.Color = [1 1 1];

% Do FDR for vS1, RSC and PPC, target:ALM
all_p_value = p_value(1,[6 7 8]);

% FDR
clear fdr_p_value_decoding
for hide = 1
    % FDR correction
    q = 0.05;
    unroll = all_p_value;
    [a,b] = sort(unroll);
    [~,r] = sort(b);
    q_threshold = r/length(all_p_value) * q;
    fdr_p_value = unroll < q_threshold;
    q_threshold = r/length(all_p_value) * 0.01;
    fdr_p_value = fdr_p_value +(unroll < q_threshold);
    q_threshold = r/length(all_p_value) * 0.001;
    fdr_p_value = fdr_p_value +(unroll < q_threshold);
end
end
