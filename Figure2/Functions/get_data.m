%% Get activity and project onto sensory and choice axes.
function data = get_data(parameters)
root_folder = [parameters.root_folder,'/Raw'] ; % CD to directory containing dataset

for animal_num = 1:7
    clearvars -except animal_list animal_num base_folder parameters root_folder data
    
    % Get animal number
    animal_struct_name = ['animal_',num2str(animal_num)];
    bpod_data = load('behaviour_data.mat', animal_struct_name);
    wavesurfer_data = load('wavesurfer_data.mat', animal_struct_name);
    imaging_data = load('imaging_data.mat', animal_struct_name);
    
    for date_num = 1:6
        clearvars -except animal_list animal_num animal_ID date_list date_num base_folder parameters  root_folder data bpod_data wavesurfer_data imaging_data animal_struct_name data
        
        %% Intialise data
        
        % Load session specific data
        session_struct_name = ['session_',num2str(date_num)];
        SessionData = bpod_data.(animal_struct_name).(session_struct_name);
        
        % Load wavesurfer data
        stim = wavesurfer_data.(animal_struct_name).(session_struct_name).trial;
        go_cue = wavesurfer_data.(animal_struct_name).(session_struct_name).go_cue;
        left_lick = wavesurfer_data.(animal_struct_name).(session_struct_name).lick_left;
        right_lick = wavesurfer_data.(animal_struct_name).(session_struct_name).lick_right;
        img_frame = wavesurfer_data.(animal_struct_name).(session_struct_name).img_frame;
        
        % Correct trials.
        for trial_type = 1:8
            correct_trial_num{trial_type} = find(SessionData.TrialTypes == trial_type & SessionData.Outcomes == 1);
            incorrect_trial_num{trial_type} = find(SessionData.TrialTypes == trial_type & SessionData.Outcomes == 0);
        end
        
        % Sampling frequency of the DAQ.
        fs_behavior = 20000;
        
        % Sampling frequency of the scanimage.
        fs_image = 9.35211;
        
        % Adjust behavior time to imaging time.
        downsamp = fs_behavior/fs_image;
        
        % Determine trial beginning and end.
        thresh = 2.5;
        stim_str = stim > thresh; % Binarize.
        stim_begin = strfind(stim_str',[0,1]) + 1;
        stim_end = strfind(stim_str',[1,0]);
        
        % Extract trial onset.
        idx = (stim_end - stim_begin) > 10000;
        stim_begin = stim_begin(idx);
        stim_end = stim_end(idx);
        
        img_frame_str = img_frame > thresh; % Binarize.
        img_frame_begin = strfind(img_frame_str',[0,1]) + 1;
        img_frame_end = strfind(img_frame_str',[1,0]);
        if img_frame_begin(1) > 300 % In case the first imaging frame onset starts after 0 time point.
            img_frame_begin = [1,img_frame_begin + 1]; % Include the first frame and shift.
        end
        img_frame_end = img_frame_end + 1; % Shift.
        % If the last imaging frame is not finished, remove it.
        if img_frame_begin(end) > img_frame_end(end)
            img_frame_begin = img_frame_begin(1:end - 1); % Remove the last frame.
        end
        
        % Determine lick bout.
        thresh = 1;
        left_lick_str = left_lick > thresh; % Binarize.
        left_lick_begin = strfind(left_lick_str',[0,1]) + 1;
        right_lick_str = right_lick > thresh; % Binarize.
        right_lick_begin = strfind(right_lick_str',[0,1]) + 1;
        
        % Down-sample behavior parameters by taking the mean and match with the imaging time.
        img_frame_middle = (img_frame_begin + img_frame_end)/2;
        for trial_num = 1:length(stim_begin)
            trial_begin_img_temp = find(abs(img_frame_middle - stim_begin(trial_num)) == min(abs(img_frame_middle - stim_begin(trial_num))));
            trial_begin_img(trial_num) = trial_begin_img_temp(1); % In case there are 2 of these, pick the first one.
            trial_end_img_temp = find(abs(img_frame_middle - stim_end(trial_num)) == min(abs(img_frame_middle - stim_end(trial_num))));
            trial_end_img(trial_num) = trial_end_img_temp(1);  % In case there are 2 of these, pick the first one.
            clear trial_begin_img_temp trial_end_img_temp
        end
        
        % Get trial on/off imaging frames.
        trial_on_img_all = [];
        for trial_num = 1:length(stim_begin)
            trial_on_img{trial_num} = [trial_begin_img(trial_num):trial_end_img(trial_num)];
            trial_on_img_all = [trial_on_img_all,trial_on_img{trial_num}];
        end
        all_img = [1:length(img_frame_begin)];
        trial_off_img_all = all_img(~ismember(all_img,trial_on_img_all));
        
        %% Sensory mode.
        trial_num_left = find(SessionData.TrialTypes == 1);
        trial_num_right = find(SessionData.TrialTypes == 2);
        
        % Load imaging data.
        for region_num = 1:8
            clear ops spks iscell cell_temp1 cell_temp2 cell zscored_spks
            
            % Load imaging data
            region_struct_name = ['region_',num2str(region_num)];
            output = imaging_data.(animal_struct_name).(session_struct_name).(region_struct_name);
            ops = output.ops;
            spks = output.spks;
            iscell = output.iscell;
            clear output
            
            cell_temp1 = find(iscell(:,1) == 1); % Select good cells.
            imaging_time = size(spks,2)/fs_image;
            
            for cell_num = 1:size(spks,1)
                zscored_spks(cell_num,:) = zscore(spks(cell_num,:));
            end
            
            % Criteria for including cells
            cell_temp2 = find(sum(zscored_spks > 10,2) >= 6*imaging_time/3600); % At least 6 times per hour zscored deconvolved signal goes above 10.
            cell = intersect(cell_temp1,cell_temp2);
            
            all_resp{region_num} = [];
            left_sensory_resp{region_num} = [];
            right_sensory_resp{region_num} = [];
            trial_averaged_all_resp{region_num} = [];
            trial_averaged_left_sensory_resp{region_num} = [];
            trial_averaged_right_sensory_resp{region_num} = [];
            for cell_num = 1:length(cell)
                for trial_num = 1:SessionData.nTrials
                    all_resp{region_num}(cell_num,trial_num,:) = smooth(zscored_spks(cell(cell_num),(trial_begin_img(trial_num) - round(fs_image*1)):(trial_begin_img(trial_num) + round(fs_image*5))),fs_image./2);
                end
                for trial_num = 1:length(trial_num_left)
                    left_sensory_resp{region_num}(cell_num,trial_num,:) = smooth(zscored_spks(cell(cell_num),(trial_begin_img(trial_num_left(trial_num)) - round(fs_image*1)):(trial_begin_img(trial_num_left(trial_num)) + round(fs_image*5))),fs_image./2);
                end
                for trial_num = 1:length(trial_num_right)
                    right_sensory_resp{region_num}(cell_num,trial_num,:) = smooth(zscored_spks(cell(cell_num),(trial_begin_img(trial_num_right(trial_num)) - round(fs_image*1)):(trial_begin_img(trial_num_right(trial_num)) + round(fs_image*5))),fs_image./2);
                end
                
                if ~isempty(all_resp{region_num}) == 1
                    trial_averaged_all_resp{region_num}(cell_num,:) = squeeze(mean(all_resp{region_num}(cell_num,:,:),2));
                else
                    trial_averaged_all_resp{region_num} = nan(length(cell),round(fs_image*6));
                end
                if ~isempty(left_sensory_resp{region_num}) == 1
                    trial_averaged_left_sensory_resp{region_num}(cell_num,:) = squeeze(mean(left_sensory_resp{region_num}(cell_num,:,:),2));
                else
                    trial_averaged_left_sensory_resp{region_num} = nan(length(cell),round(fs_image*6));
                end
                if ~isempty(right_sensory_resp{region_num}) == 1
                    trial_averaged_right_sensory_resp{region_num}(cell_num,:) = squeeze(mean(right_sensory_resp{region_num}(cell_num,:,:),2));
                else
                    trial_averaged_right_sensory_resp{region_num} = nan(length(cell),round(fs_image*6));
                end
            end
        end
        
        % Sensory mode
        for region_num = 1:8
            if size(trial_averaged_all_resp{region_num},1) >= 20 % At least 20 cells in a given region.
                sensory_diff{region_num} = mean(trial_averaged_right_sensory_resp{region_num}(:,round(fs_image*1):round(fs_image*2)),2) - mean(trial_averaged_left_sensory_resp{region_num}(:,round(fs_image*1):round(fs_image*2)),2);
                norm_sensory_diff{region_num} = sensory_diff{region_num}./((sum(sensory_diff{region_num}.^2)).^0.5);
                
                % Project to sensory mode.
                for trial_num = 1:SessionData.nTrials
                    sensory_mode{region_num}(trial_num,:) = (squeeze(all_resp{region_num}(:,trial_num,:)))'*norm_sensory_diff{region_num};
                end
            else
                sensory_mode{region_num} = nan(SessionData.nTrials,round(fs_image*6));
            end
        end
        
        %% Choice mode
        choice_left = find(SessionData.TrialTypes == 1 & SessionData.Outcomes == 1 | SessionData.TrialTypes == 2 & SessionData.Outcomes == 0);
        choice_right = find(SessionData.TrialTypes == 1 & SessionData.Outcomes == 0 | SessionData.TrialTypes == 2 & SessionData.Outcomes == 1);
        
        % Load imaging data.
        for region_num = 1:8
            clear ops spks iscell cell_temp1 cell_temp2 cell zscored_spks
            
            % Load imaging data
            region_struct_name = ['region_',num2str(region_num)];
            output = imaging_data.(animal_struct_name).(session_struct_name).(region_struct_name);
            ops = output.ops;
            spks = output.spks;
            iscell = output.iscell;
            clear output
            
            cell_temp1 = find(iscell(:,1) == 1); % Select good cells.
            imaging_time = size(spks,2)/fs_image;
            
            for cell_num = 1:size(spks,1)
                zscored_spks(cell_num,:) = zscore(spks(cell_num,:));
            end
            cell_temp2 = find(sum(zscored_spks > 10,2) >= 6*imaging_time/3600); % At least 6 times per hour zscored deconvolved signal goes above 10.
            cell = intersect(cell_temp1,cell_temp2);
            
            all_resp{region_num} = [];
            left_choice_resp{region_num} = [];
            right_choice_resp{region_num} = [];
            trial_averaged_all_resp{region_num} = [];
            trial_averaged_left_choice_resp{region_num} = [];
            trial_averaged_right_choice_resp{region_num} = [];
            for cell_num = 1:length(cell)
                for trial_num = 1:SessionData.nTrials
                    all_resp{region_num}(cell_num,trial_num,:) = smooth(zscored_spks(cell(cell_num),(trial_begin_img(trial_num) - round(fs_image*1)):(trial_begin_img(trial_num) + round(fs_image*5))),fs_image./2);
                end
                for trial_num = 1:length(choice_left)
                    left_choice_resp{region_num}(cell_num,trial_num,:) = smooth(zscored_spks(cell(cell_num),(trial_begin_img(choice_left(trial_num)) - round(fs_image*1)):(trial_begin_img(choice_left(trial_num)) + round(fs_image*5))),fs_image./2);
                end
                for trial_num = 1:length(choice_right)
                    % Cell x trial x frames
                    right_choice_resp{region_num}(cell_num,trial_num,:) = smooth(zscored_spks(cell(cell_num),(trial_begin_img(choice_right(trial_num)) - round(fs_image*1)):(trial_begin_img(choice_right(trial_num)) + round(fs_image*5))),fs_image./2);
                end
                
                if ~isempty(all_resp{region_num}) == 1
                    trial_averaged_all_resp{region_num}(cell_num,:) = squeeze(mean(all_resp{region_num}(cell_num,:,:),2));
                else
                    trial_averaged_all_resp{region_num} = nan(length(cell),round(fs_image*6));
                end
                if ~isempty(left_choice_resp{region_num}) == 1
                    trial_averaged_left_choice_resp{region_num}(cell_num,:) = squeeze(mean(left_choice_resp{region_num}(cell_num,:,:),2));
                else
                    trial_averaged_left_choice_resp{region_num} = nan(length(cell),round(fs_image*6));
                end
                if ~isempty(right_choice_resp{region_num}) == 1
                    %Cell x frame
                    trial_averaged_right_choice_resp{region_num}(cell_num,:) = squeeze(mean(right_choice_resp{region_num}(cell_num,:,:),2));
                else
                    trial_averaged_right_choice_resp{region_num} = nan(length(cell),round(fs_image*6));
                end
            end
        end
        
        for region_num = 1:8
            if size(trial_averaged_all_resp{region_num},1) >= 20 % At least 20 cells in a given region.
                choice_diff{region_num} = mean(trial_averaged_right_choice_resp{region_num}(:,round(fs_image*3.5):round(fs_image*4)),2) - mean(trial_averaged_left_choice_resp{region_num}(:,round(fs_image*3.5):round(fs_image*4)),2);
                norm_choice_diff{region_num} = choice_diff{region_num}./((sum(choice_diff{region_num}.^2)).^0.5);
                
                % Project to choice mode.
                for trial_num = 1:SessionData.nTrials
                    choice_mode{region_num}(trial_num,:) = (squeeze(all_resp{region_num}(:,trial_num,:)))'*norm_choice_diff{region_num};
                end
            else
                choice_mode{region_num} = nan(SessionData.nTrials,round(fs_image*6));
            end
        end
        
        %% Gram-Schmidt process to orthogonalize these modes on a trial-by-trial basis.
        for region_num = 1:8
            if ~isempty(choice_mode{region_num}) == 1
                for trial_num = 1:SessionData.nTrials
                    clear X d n m R Q D
                    X = [sensory_mode{region_num}(trial_num,:);choice_mode{region_num}(trial_num,:)]';
                    
                    %%%%% Written by Mo Chen (sth4nth@gmail.com). %%%%%
                    [d,n] = size(X);
                    m = min(d,n);
                    R = eye(m,n);
                    Q = zeros(d,m);
                    D = zeros(1,m);
                    for i = 1:m
                        R(1:i-1,i) = bsxfun(@times,Q(:,1:i-1),1./D(1:i-1))'*X(:,i);
                        Q(:,i) = X(:,i)-Q(:,1:i-1)*R(1:i-1,i);
                        D(i) = dot(Q(:,i),Q(:,i));
                    end
                    R(:,m+1:n) = bsxfun(@times,Q,1./D)'*X(:,m+1:n);
                    %%%%% Written by Mo Chen (sth4nth@gmail.com). %%%%%
                    
                    Q_trial{region_num}{trial_num} = Q;
                end
                
                choice_mode_GS{region_num} = [];
                sensory_mode_GS{region_num} = [];
                for trial_num = 1:SessionData.nTrials
                    choice_mode_GS{region_num}(trial_num,:) = Q_trial{region_num}{trial_num}(:,2);
                    sensory_mode_GS{region_num}(trial_num,:) = Q_trial{region_num}{trial_num}(:,1);
                end
            else
                choice_mode_GS{region_num} = nan(SessionData.nTrials,round(fs_image*6));
                sensory_mode_GS{region_num} = nan(SessionData.nTrials,round(fs_image*6));
            end
        end
        
        %% Save processed data into data_field
        data.all_choice_mode_GS{animal_num}{date_num} = choice_mode_GS;
        data.all_sensory_mode_GS{animal_num}{date_num} = sensory_mode_GS;
        data.all_trial_averaged_all_resp{animal_num}{date_num} = trial_averaged_all_resp;
        data.all_trial_averaged_left_choice_resp{animal_num}{date_num} = trial_averaged_left_choice_resp;
        data.all_trial_averaged_right_choice_resp{animal_num}{date_num} =trial_averaged_right_choice_resp;
        data.all_trial_averaged_left_sensory_resp{animal_num}{date_num} = trial_averaged_left_sensory_resp;
        data.all_trial_averaged_right_sensory_resp{animal_num}{date_num} = trial_averaged_right_sensory_resp;
        data.all_left_choice_resp{animal_num}{date_num} = left_choice_resp;
        data.all_right_choice_resp{animal_num}{date_num} = right_choice_resp;
        data.all_left_sensory_resp{animal_num}{date_num} = left_sensory_resp;
        data.all_right_sensory_resp{animal_num}{date_num} = right_sensory_resp;
        
        data.all_SessionData{animal_num}{date_num} = SessionData;
        
        data.choice_left{animal_num}{date_num} = choice_left;
        data.choice_right{animal_num}{date_num} = choice_right;
        data.trial_num_left{animal_num}{date_num} = trial_num_left;
        data.trial_num_right{animal_num}{date_num} = trial_num_right;
        data.all_trial_type{animal_num}{date_num} = SessionData.TrialTypes;
        
        data.norm_choice_diff{animal_num}{date_num} = norm_choice_diff;
        data.norm_sensory_diff{animal_num}{date_num} = norm_sensory_diff;
        
    end
end
end