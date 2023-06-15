%% To plot choice mode main figure for 1 example trained RNN
function RNN_data = get_example_choice_mode_RNN(folder_path,matfile_folder_name)
cd(folder_path)
num_units = 807;
file_list = dir;
file_list = file_list(contains({file_list.name},'expert')); % Extract files containing animal ID.

% Loop over different RNNs
for k = 1 % k = 1 demo session
    seed_file_list = file_list([1 2]);
    for i = 1:size(seed_file_list,1)
        % Trial activity
        fs_image = 93.5;
        pre_resp_epoch = [1:round(fs_image*3.5)];
        delay_epoch = [round(fs_image*2.5):round(fs_image*3.5)];
        post_sample_epoch = [round(fs_image*1.5):round(fs_image*3.5)];
        choice_epoch = [round(fs_image*3.0):round(fs_image*3.5)];
        sample_epoch = round(fs_image*0.5) : round(fs_image*1.5);
        ITI_epoch = 1:round(fs_image*0.5);
        
        %% Get MSE of each seed using target and test activity
        cd([folder_path,'/',seed_file_list(i).name,]);
        load readouts_r.mat
        fit_R = R;
        fit_L = L;
        clear R L
        
        load ground_truth_for_r.mat
        target_R = R;
        target_L = L;
        clear R L
        
        RNN_data.MSE_right(k,i) = mean(mean((target_R - fit_R).^2));
        RNN_data.MSE_left(k,i) = mean(mean((target_L - fit_L).^2));

        %% Get choice mode from perturbed activity
        cd([folder_path,'/',seed_file_list(i).name,'/',matfile_folder_name]);
        file_name = dir;
        file_name = file_name(contains({file_name.name},'trials')); % Extract files containing animal ID.
        
        load(file_name.name);
        
        clear left_choice_mode_norm right_choice_mode_norm L_pos_choice_mode_norm R_neg_choice_mode_norm
        
        averaged_R = squeeze(mean(R,1));
        averaged_L = squeeze(mean(L,1));
        
        clear norm_choice_diff choice_diff
        choice_diff = squeeze(mean(averaged_R(:,choice_epoch),2) - mean(averaged_L(:,choice_epoch),2));
        
        norm_choice_diff = choice_diff./((sum(choice_diff.^2)).^0.5);
        norm_choice_diff = norm_choice_diff(1:num_units);
        
        for trial_num = 1:size(R,1)
            trial_R = squeeze(R(trial_num,1:num_units,:));
            trial_L = squeeze(L(trial_num,1:num_units,:));
            trial_L_pos = squeeze(L_pos(trial_num,1:num_units,:));
            
            % Set criteria to remove abberant trials with std > 6x median STD
            if (sum(mean(trial_R(:,delay_epoch),1)> (6*median(std(trial_R,1)))) < 1) & ( sum( mean(trial_L(:,delay_epoch),1) > (6*mean(std(trial_L,1))) ) < 1)
                
                % Project to sensory mode.
                left_choice_mode = trial_L'*norm_choice_diff;
                right_choice_mode = trial_R'*norm_choice_diff;
                L_pos_choice_mode = trial_L_pos'*norm_choice_diff;
                
                % Substract baseline
                left_choice_mode_substracted = left_choice_mode - nanmean(left_choice_mode(ITI_epoch));
                right_choice_mode_substracted = right_choice_mode- nanmean(right_choice_mode(ITI_epoch));
                L_pos_choice_mode_substracted = L_pos_choice_mode - nanmean(L_pos_choice_mode(ITI_epoch));
                
                if max(right_choice_mode_substracted(pre_resp_epoch)) > max(abs(left_choice_mode_substracted(pre_resp_epoch)))
                    max_value = max(right_choice_mode_substracted(pre_resp_epoch));
                elseif max(right_choice_mode_substracted(pre_resp_epoch)) < max(abs(left_choice_mode_substracted(pre_resp_epoch)))
                    max_value = abs(min(left_choice_mode_substracted(pre_resp_epoch)));
                end
                
                left_choice_mode_std_norm(trial_num,:) = left_choice_mode_substracted./ max(right_choice_mode_substracted(pre_resp_epoch));
                right_choice_mode_std_norm(trial_num,:) = right_choice_mode_substracted./ max(right_choice_mode_substracted(pre_resp_epoch));
                L_pos_choice_mode_std_norm(trial_num,:) = L_pos_choice_mode_substracted./ max(right_choice_mode_substracted(pre_resp_epoch));
                
                if (L_pos_choice_mode_std_norm(trial_num,end) > left_choice_mode_std_norm(trial_num,327) )  & (right_choice_mode_std_norm(trial_num,327) > left_choice_mode_std_norm(trial_num,327))
                    L_pos_mid_point(trial_num) = (L_pos_choice_mode_std_norm(trial_num,327) - left_choice_mode_std_norm(trial_num,327) ) / (right_choice_mode_std_norm(trial_num,327) - left_choice_mode_std_norm(trial_num,327));
                else
                    L_pos_mid_point(trial_num) = 0;
                end
                
                if L_pos_mid_point(trial_num) > 1
                    L_pos_mid_point(trial_num) = 1;
                end
            else
                L_pos_mid_point(trial_num) = nan;
            end
        end
        
        if i == 1
            RNN_data.PPC.choice_mode.left_choice_mode_std_norm{k} = left_choice_mode_std_norm;
            RNN_data.PPC.choice_mode.right_choice_mode_std_norm{k} = right_choice_mode_std_norm;
            RNN_data.PPC.choice_mode.L_pos_choice_mode_std_norm{k} = L_pos_choice_mode_std_norm;
            RNN_data.PPC.choice_mode.L_pos_mid_point{k} = L_pos_mid_point;

        elseif i == 2
            RNN_data.vS1.choice_mode.left_choice_mode_std_norm{k} = left_choice_mode_std_norm;
            RNN_data.vS1.choice_mode.right_choice_mode_std_norm{k} = right_choice_mode_std_norm;
            RNN_data.vS1.choice_mode.L_pos_choice_mode_std_norm{k} = L_pos_choice_mode_std_norm;
            RNN_data.vS1.choice_mode.L_pos_mid_point{k} = L_pos_mid_point;
        end
        all_proportion_switch_left(k,i) = sum(L_pos_mid_point>0.5);
    end
end
RNN_data.all_proportion_switch_left = all_proportion_switch_left;
end