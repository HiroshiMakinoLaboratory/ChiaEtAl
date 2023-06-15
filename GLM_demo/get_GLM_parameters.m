%% (Demo) Generalized linear model (GLM).
% For GLM cross-validation, cvpartitionImpl.m in /Applications/MATLAB_R2018b.app/toolbox/stats/stats/+internal/+stats line 113 needs to be changed to "properties(GetAccess = public, SetAccess = public)".
parameters.root_folder = uigetdir;
cd([parameters.root_folder,'/GLM_demo'])

load example_bpod_data.mat
load example_roimatchpub.mat
load example_suite2p_data.mat
load example_wavesurfer_data.mat

fs_image = 9.352;
x_size = 500;
y_size = 500;
region = 8;
xy_pixel_ratio = 0.4;
suite2p_param_aspect = 2.400;
suite2p_param_diameter = [12 5];
suite2p_param_threshold_scaling = 0.5000;
region_num = 8; % For example region
session_num = 3; % For example expert session

total_num_cells_session = 0;

% Get the order of the sessions from roiMatchPub file (obtained from the path of the image) (In case the sessions are not in order)
for k = 1:length(roiMatchData.allRois)
    idx = strfind(lower(roiMatchData.allRois{k}),'imaging');
    session_name(k,:) = str2num(roiMatchData.allRois{k}(idx+8 :idx+13));
    clear idx
end

% Additional step to remove unreliable matched cells and/or duplicates from roiMatchPub file
% Valid Match: Ensure that at least more than 2 sessions matches i.e.
% e.g. Sess1: 44-22-28 match, Sess2: 44-22-38 match and Sess3: 44-22-38 match. (match--> above threshold)
% Invalid Match : Sess1: 44-22-NA match, Sess2: 44-22-38 match and Sess3: 44-22-38 match
allsession_matched = roiMatchData.allSessionMapping;
allsession_matched(:,end+1) = zeros(size(allsession_matched,1),1);
usecol = 1; % Using first column/session as a reference
number_sessions = size(roiMatchData.allSessionMapping,2);

for irow = 1:length(roiMatchData.allSessionMapping(:,usecol))
    % Excluding the first row
    % In roiMatchData.mapping, the pairwise comparisons are matched repeatedly for all sessions (i.e. if there are 3 sessions, the number of matches will be repeated row-wise.
    % Hence, if you can find the same matching 3 times, it means that the same cell is matched all 3 times for all sessions)
    match3idx = find( sum(roiMatchData.mapping(:,setdiff(1:number_sessions,usecol)) == roiMatchData.allSessionMapping(irow,setdiff(1:number_sessions,usecol)),2) == size(roiMatchData.allSessionMapping,2)-1);
    if length(match3idx) >= (size(roiMatchData.allSessionMapping,2))
        allsession_matched(irow,end) = 1; % If it is 1, it means that the matching is valid where there all 3 sessions the same cells were matched.
    end
    clear match3idx
end

% Get idx of all cells that are matched with the other sessions
match_cell_index{region_num} = roiMatchData.allSessionMapping(:,session_num); % Using specific idx of specific session
% Remove non-valid matches
match_cell_index{region_num} = match_cell_index{region_num}(logical(allsession_matched(:,end)));
clear allsession_matched

% Ensure no repeated matched cells
assert(length(unique(match_cell_index{region_num})) == length(match_cell_index{region_num}),'Repeated index in matched cell index (from RoiMatchPub)')

% Get total number of cells in a session to check if number of predictors are correct
total_num_cells_session = total_num_cells_session + size(match_cell_index{region_num},1);

region_labels = [];
%clearvars  -except animal_ID date fs_image x_size y_size region xy_pixel_ratio Bpod_protocol listing region_file region_num cell x_pixel y_pixel whole_img spike_all F_raw_all
whole_img{region_num} = [];
whole_imge{region_num} = [];
spike_all{region_num} = [];
F_raw_all{region_num} = [];
valid_cell{region_num} = find(iscell(:,1) == 1); % Select good cells.

% Note that matched cells idxs are idxs for valid cells (i.e. diff idx from iscell)
% Assign pixels.

% Only taking from valid cells idx
for cell_num = 1:length(valid_cell{region_num})
    x_pixel{region_num}{cell_num} = stat{valid_cell{region_num}(cell_num)}.xpix;
    y_pixel{region_num}{cell_num} = stat{valid_cell{region_num}(cell_num)}.ypix;
end
% From valid cells variable, obtained matched cells only
for i = 1: size(match_cell_index{region_num},1) % Adjust to match index
    cell_num = match_cell_index{region_num}(i);
    matched_x_pixel{region_num}{i} = x_pixel{region_num}{cell_num};
    matched_y_pixel{region_num}{i} = y_pixel{region_num}{cell_num};
end
x_pixel = matched_x_pixel;
y_pixel = matched_y_pixel;

% Activity of ALL individual neurons.
spike = spks(valid_cell{region_num},:); % Deconvolved data.
F_raw = F(valid_cell{region_num},:) - 0.7*Fneu(valid_cell{region_num},:);
whole_img{region_num} = ops.meanImg;
whole_imge{region_num} = ops.meanImgE;

if ~isnan(match_cell_index{region_num}) % Concatenate ONLY matched cells
    spike_all{region_num} = [spike_all{region_num};spike(match_cell_index{region_num},:)]; % Matched cells will be saved according to order in roiMatchPub idx.
    F_raw_all{region_num} = [F_raw_all{region_num};F_raw(match_cell_index{region_num},:)];
end
% Create a region labels to check if the number of cells concatenated in each region is consistent
region_labels = [region_labels; repmat(region_num,size(spike_all{region_num},1),1)];
clear spike F_raw
assert(size(region_labels,1) == total_num_cells_session,'Error in number of matched cells');
%% Analyze behavior.
% Locate Bpod folders and load Bpod file.

% Tactile 2AFC DAQ setup.
trial_ch = 1;
go_cue_ch = 2;
lick_left_ch = 5;
lick_right_ch = 6;
img_frame_ch = 8;

% Read DAQ data from waversurfer.
data_fieldnames = fieldnames(data);
trial = data.(data_fieldnames{2}).analogScans(:,trial_ch);
go_cue = data.(data_fieldnames{2}).analogScans(:,go_cue_ch);
lick_left = data.(data_fieldnames{2}).analogScans(:,lick_left_ch);
lick_right = data.(data_fieldnames{2}).analogScans(:,lick_right_ch);
img_frame = data.(data_fieldnames{2}).analogScans(:,img_frame_ch);

% Get trial types from bpod data
correct_control_B(:,1) = SessionData.TrialTypes; % Will contain information on trial type and correct/incorrect trial
for trial_num_bpod = 1:SessionData.nTrials
    % Differs from pre_training A because it does not require
    if (SessionData.TrialTypes(trial_num_bpod) >= 3)
        correct_control_B(trial_num_bpod,2) = 3 ; % catch trial == 3
    elseif (~isnan(SessionData.RawEvents.Trial{trial_num_bpod}.States.Reward(1)) && SessionData.TrialTypes(trial_num_bpod) <= 2 )
        correct_control_B(trial_num_bpod,2) = 1 ; % Reward trial == 1
    elseif (~isnan(SessionData.RawEvents.Trial{trial_num_bpod}.States.Punish(1)) && (SessionData.RawEvents.Trial{trial_num_bpod}.States.Response(2) >= SessionData.TrialSettings(1).GUI.StimDuration + SessionData.TrialSettings(1).GUI.DelayDuration + SessionData.TrialSettings(1).GUI.ResponseDuration  + SessionData.TrialSettings(1).GUI.GoCueDuration) && SessionData.TrialTypes(trial_num_bpod) <= 2)
        correct_control_B(trial_num_bpod,2) = 2 ; % no response trial == 2
    end
end

% Get trial types
trial_type = nan(1,length(correct_control_B(:,1)))';
trial_type(find(correct_control_B(:,1) == 1 & correct_control_B(:,2) == 1)) = 1; % Correct left
trial_type(find(correct_control_B(:,1) == 1 & correct_control_B(:,2) == 0)) = 2;% Incorrect left
trial_type(find(correct_control_B(:,1) == 2 & correct_control_B(:,2) == 1)) = 3; % Correct right
trial_type(find(correct_control_B(:,1) == 2 & correct_control_B(:,2) == 0)) = 4; % Inorrect right
trial_type(find(correct_control_B(:,2) == 2)) =  0; % No response

% Sampling frequency of the DAQ.
fs_behavior = data.header.AcquisitionSampleRate;

% Adjust behavior time to imaging time.
downsamp = fs_behavior/fs_image;

% Determine trial begining and end.
thresh = 2.5;
trial_str = trial > thresh; % Binarize.
trial_begin = strfind(trial_str',[0,1]) + 1;
% Trial offset was set to be at the beginning of ITI.
for trial_num = 1:length(trial_begin)
    trial_end(trial_num) = trial_begin(trial_num) + fs_behavior.*SessionData.RawEvents.Trial{trial_num}.States.ITI(1); % To start of ITI
end
go_cue_str = go_cue > thresh; % Binarize.
go_cue_begin = strfind(go_cue_str',[0,1]) + 1; % Go Cue begin is the end of delay epoch
go_cue_end = strfind(go_cue_str',[1,0]);

img_frame_str = img_frame > thresh; % Binarize.
img_frame_begin = strfind(img_frame_str',[0,1]) + 1;
img_frame_end = strfind(img_frame_str',[1,0]);

% Edited - add a stimulus middle (i.e. 0.5s from trial onset)
stim_middle = trial_begin + (fs_behavior * SessionData.TrialSettings(1).GUI.StimDuration *0.5);

delay_begin = trial_begin + fs_behavior * SessionData.TrialSettings(1).GUI.StimDuration + 1;
delay_end = trial_begin + fs_behavior * SessionData.TrialSettings(1).GUI.StimDuration + fs_behavior * SessionData.TrialSettings(1).GUI.DelayDuration;
assert(SessionData.TrialSettings(1).GUI.DelayDuration == 2);
delay_middle = trial_begin + (fs_behavior * SessionData.TrialSettings(1).GUI.StimDuration) +  (fs_behavior * 0.75); % Changed to 0.75

if img_frame_begin(1) > 0.5.*(fs_behavior./fs_image) % In case the first imaging frame onset starts after 0 time point.
    img_frame_begin = [1,img_frame_begin + 1]; % Include the first frame and shift.
end
img_frame_end = img_frame_end + 1; % Shift.
if img_frame_begin(end) > img_frame_end(end)
    img_frame_begin = img_frame_begin(1:end - 1); % Remove the last frame.
end

% Down-sample behavior parameters by taking the mean and match with the imaging time.
img_frame_middle = (img_frame_begin + img_frame_end)/2;
for trial_num = 1:length(trial_begin)
    
    trial_begin_img_temp = find(abs(img_frame_middle - trial_begin(trial_num)) == min(abs(img_frame_middle - trial_begin(trial_num))));
    trial_begin_img(trial_num) = trial_begin_img_temp(1); % In case there are 2 of these, pick the first one.
    trial_end_img_temp = find(abs(img_frame_middle - trial_end(trial_num)) == min(abs(img_frame_middle - trial_end(trial_num))));
    trial_end_img(trial_num) = trial_end_img_temp(1); % In case there are 2 of these, pick the first one.
    
    % Edited to add tact_stim middle
    stim_middle_img_temp = find(abs(img_frame_middle - stim_middle(trial_num)) == min(abs(img_frame_middle - stim_middle(trial_num))));
    stim_middle_img(trial_num) = stim_middle_img_temp(1); % In case there are 2 of these, pick the first one.
    
    delay_begin_img_temp = find(abs(img_frame_middle - delay_begin(trial_num)) == min(abs(img_frame_middle - delay_begin(trial_num))));
    delay_begin_img(trial_num) = delay_begin_img_temp(1); % In case there are 2 of these, pick the first one.
    delay_end_img_temp = find(abs(img_frame_middle - delay_end(trial_num)) == min(abs(img_frame_middle - delay_end(trial_num))));
    delay_end_img(trial_num) = delay_end_img_temp(1); % In case there are 2 of these, pick the first one.
    delay_middle_img_temp = find(abs(img_frame_middle - delay_middle(trial_num)) == min(abs(img_frame_middle - delay_middle(trial_num))));
    delay_middle_img(trial_num) = delay_middle_img_temp(1); % In case there are 2 of these, pick the first one.
    
    go_cue_begin_img_temp = find(abs(img_frame_middle - go_cue_begin(trial_num)) == min(abs(img_frame_middle - go_cue_begin(trial_num))));
    go_cue_begin_img(trial_num) = go_cue_begin_img_temp(1); % In case there are 2 of these, pick the first one.
    
    clear trial_begin_img_temp trial_end_img_temp delay_begin_img_temp delay_end_img_temp delay_middle_img_temp
end

% Determine lick bout.
% Adjust threshold according to animal (e.g. XW0077-L keeps touch lickport - alot of noise )
thresh = 2.5;
lick_left_str = lick_left > thresh; % Binarize.
lick_right_str = lick_right > thresh; % Binarize.
lick_left_begin = strfind(lick_left_str',[0,1]) + 1;
lick_right_begin = strfind(lick_right_str',[0,1]) + 1;

% Get the first lick bout after the go cue end.
for trial_num = 1:length(trial_begin)
    %trial_lick{trial_num} = trial_begin(trial_num) < lick_center_begin & lick_center_begin < trial_end(trial_num);
    trial_lick_left{trial_num} = trial_begin(trial_num) < lick_left_begin & lick_left_begin < trial_end(trial_num); % Find all the lick bout that happens between trial start and end
    trial_lick_right{trial_num} = trial_begin(trial_num) < lick_right_begin & lick_right_begin < trial_end(trial_num);
    %
    
    % Find the first lick bout for left lick within trial epoch
    if ~isempty(min(strfind(trial_lick_left{trial_num},[1]))) == 1 %
        %post_go_cue_lick_center_idx(trial_num) = min(strfind(trial_lick{trial_num},1));
        post_go_cue_lick_left_idx(trial_num) = min(strfind(trial_lick_left{trial_num},1));
    else
        %post_go_cue_lick_center_idx(trial_num) = nan;
        post_go_cue_lick_left_idx(trial_num) = nan;
    end
    
    if ~isempty(min(strfind(trial_lick_right{trial_num},[1]))) == 1 %
        %post_go_cue_lick_center_idx(trial_num) = min(strfind(trial_lick{trial_num},1));
        post_go_cue_lick_right_idx(trial_num) = min(strfind(trial_lick_right{trial_num},1));
    else
        %post_go_cue_lick_center_idx(trial_num) = nan;
        post_go_cue_lick_right_idx(trial_num) = nan;
    end
end

% post_go_cue_lick_right_idx is the frame number idx for the first left/right lick bout in WS framerate
% Assign 1 for lick event to nearby imaging frames.
for trial_num = 1:length(trial_begin)
    if isnan(post_go_cue_lick_left_idx(trial_num)) == 1
        post_go_cue_lick_left_img(trial_num) = nan;
    else
        % Added a min to min(find(abs(img...))) because it was possible to get >1 lick event to nearby imaging frame
        post_go_cue_lick_left_img(trial_num) = min(find(abs(img_frame_middle - lick_left_begin(post_go_cue_lick_left_idx(trial_num))) == min(abs(img_frame_middle - lick_left_begin(post_go_cue_lick_left_idx(trial_num))))));
    end
    
    if isnan(post_go_cue_lick_right_idx(trial_num)) == 1
        post_go_cue_lick_right_img(trial_num) = nan;
    else
        % Added a min to min(find(abs(img...))) because it was possible to get >1 lick event to nearby imaging frame
        post_go_cue_lick_right_img(trial_num) = min(find(abs(img_frame_middle - lick_right_begin(post_go_cue_lick_right_idx(trial_num))) == min(abs(img_frame_middle - lick_right_begin(post_go_cue_lick_right_idx(trial_num))))));
    end
    
    assert(trial_begin_img(trial_num)<post_go_cue_lick_right_img(trial_num)<trial_end_img(trial_num))
    assert(trial_begin_img(trial_num)<post_go_cue_lick_left_img(trial_num)<trial_end_img(trial_num))
    
    % Run assertion to ensure that it is correct for right trials
    if ~isnan(post_go_cue_lick_right_img(trial_num)) & ~isnan(post_go_cue_lick_left_img(trial_num))
        if trial_type(trial_num)==3
            if ~((post_go_cue_lick_right_img(trial_num)<post_go_cue_lick_left_img(trial_num)))
                warning('lick trial onset error')
            end
        end
    end
    % Run assertion to ensure that it is correct for left trials
    if ~isnan(post_go_cue_lick_right_img(trial_num)) & ~isnan(post_go_cue_lick_left_img(trial_num))
        if trial_type(trial_num)==1
            if ~((post_go_cue_lick_right_img(trial_num)>post_go_cue_lick_left_img(trial_num)))
                warning('lick trial onset error')
            end
        end
    end
end

%% Analyse movement using .csv file from deeplabcut, note: movie frame rate has to be downsampled
cd([parameters.root_folder,'/GLM_demo/DLC'])
listing = dir;
csv_file_list = listing(contains({listing.name},'.csv')); % Extract files containing .csv
csv_file_list = csv_file_list(~contains({csv_file_list.name},'._')); % For some reason, duplicate hidden files with ._ can be created unknowingly

% Start sorting of all csv file names
name_position_vector = [];
for i = 1: length(csv_file_list)
    start_position = strfind(csv_file_list(i).name,'DLC');
    index = str2num(csv_file_list(i).name(start_position-4:start_position-1));
    name_position_vector = [name_position_vector; index];
end

name_position_vector = name_position_vector +1; % To account for indexing starting at zero
[~, idx] = sortrows((name_position_vector));
sorted_csv_file_list = csv_file_list(idx); % Sort files based on the name_position_vector index

% Load vector for trial onset
load TrialOnsetVector.mat

% Run quick checks to ensure that the files are correct
% Check the trial onset vector obtained from running gettrialonset.mat
if length(find(trialonset_vector_adjusted == 1)) ~= SessionData.nTrials % Display warning of number of trials in movie differs from the supposed number of trials in a session (i.e. 180)
    disp('ERROR - Number of trials detected from movie differs from number of trials')
end

if length(name_position_vector) ~= max(name_position_vector)
    disp('ERROR - Number of DLC .csv files do not match the number of videos')
end

% Combine all the CSV files within the directory to form a combined table
% Concatenate all the CSV files together to form a vector
csv_table = [];
combined_csv_table=[];
for csv_num = 1:length(sorted_csv_file_list)
    csv_table = readtable(sorted_csv_file_list(csv_num).name, 'HeaderLines',3);
    combined_csv_table = [combined_csv_table; csv_table];
end

% Convert deeplabcut output CSV table into matrix
combined_csv_table = combined_csv_table{:,:};
% Add column trialonset_vector_adjusted to the last column
combined_csv_table = [combined_csv_table trialonset_vector_adjusted'];

% Get trial onset and offset for movie
movie_trial_begin = find(combined_csv_table(:,end)==1)';
% Trial offset was set to be at the beginning of ITI.
for trial_num = 1:length(movie_trial_begin)
    movie_trial_end(trial_num) = movie_trial_begin(trial_num) + round(video_frame_rate).*SessionData.RawEvents.Trial{trial_num}.States.ITI(1); % Trial ends at start of ITI
end

% Extract selected forepaw movement from matrix
rightforepaw_session_movement = combined_csv_table(:,2:3)';
leftforepaw_session_movement = combined_csv_table(:,5:6)';

% likelihood of position from matrix
rightforepaw_session_likelihood = combined_csv_table(:,4)';
leftforepaw_session_likelihood = combined_csv_table(:,7)';

% Processing x and y coordinates data from DLC
% First step: Removing occluded or low confidence points
% If Likelihood of deeplabcut analysis is < 0.95 use the mean of previous and next frame with >0.95 to replace
right_bad_frames = find(rightforepaw_session_likelihood<0.95);
right_good_frames = find(rightforepaw_session_likelihood>=0.95);
for i = 1: length(right_bad_frames)
    bad_frame = right_bad_frames(i);
    neg_diff = (right_good_frames - bad_frame);
    pos_only_neg_diff = neg_diff(find(neg_diff>0));
    [minimum, ~] = min(pos_only_neg_diff(1,:));
    neg_closestvalue = bad_frame + minimum;
    clear minimum
    
    pos_diff = (bad_frame - right_good_frames.');
    pos_only_pos_diff = pos_diff(find(pos_diff>0));
    [minimum, ~] = min(pos_only_pos_diff);
    pos_closestvalue = bad_frame - minimum;
    
    % If either ends are unavailable, i.e. all the frames are <0.95, use only the closest frame on the other end
    if isempty(neg_closestvalue)
        neg_closestvalue =  pos_closestvalue;
        % Check if code works
        assert(rightforepaw_session_likelihood(pos_closestvalue) >= 0.95,'Error with DLC point processing')
        if pos_closestvalue+1 < length(rightforepaw_session_likelihood)
            assert(rightforepaw_session_likelihood(pos_closestvalue+1) < 0.95,'Error with DLC point processing')
        end
    end
    if isempty(pos_closestvalue)
        pos_closestvalue =  neg_closestvalue;
        
        % Check if code works
        assert(rightforepaw_session_likelihood(neg_closestvalue) >= 0.95,'Error with DLC point processing')
        if neg_closestvalue-1 > 0
            assert(rightforepaw_session_likelihood(neg_closestvalue-1) < 0.95,'Error with DLC point processing')
        end
    end
    
    % Replace the bad outlier frame x and y values with mean of closest frame values
    rightforepaw_session_movement(:,bad_frame) = (rightforepaw_session_movement(:,neg_closestvalue) + rightforepaw_session_movement(:,pos_closestvalue))/2;
    
    clear neg_diff pos_only_neg_diff neg_closestvalue pos_diff pos_only_npos_diff pos_closestvalue
end

% If Likelihood of deeplabcut analysis is < 0.95 use the mean of previous and next frame with >0.95 to replace
left_bad_frames = find(leftforepaw_session_likelihood<0.95);
left_good_frames = find(leftforepaw_session_likelihood>=0.95);
for i = 1: length(left_bad_frames)
    bad_frame = left_bad_frames(i);
    neg_diff = (left_good_frames - bad_frame);
    pos_only_neg_diff = neg_diff(find(neg_diff>0));
    [minimum, ~] = min(pos_only_neg_diff(1,:));
    neg_closestvalue = bad_frame + minimum;
    clear minimum
    
    pos_diff = (bad_frame - left_good_frames.');
    pos_only_pos_diff = pos_diff(find(pos_diff>0));
    [minimum, ~] = min(pos_only_pos_diff);
    pos_closestvalue = bad_frame - minimum;
    
    % If either ends are unavailable, i.e. all the frames are <0.95, use only the closest frame on the other end
    if isempty(neg_closestvalue)
        neg_closestvalue =  pos_closestvalue;
        % Check if code works
        assert(leftforepaw_session_likelihood(pos_closestvalue) >= 0.95,'Error with DLC point processing')
        if pos_closestvalue+1 < length(leftforepaw_session_likelihood)
            assert(leftforepaw_session_likelihood(pos_closestvalue+1) < 0.95,'Error with DLC point processing')
        end
    end
    if isempty(pos_closestvalue)
        pos_closestvalue =  neg_closestvalue;
        
        % Check if code works
        assert(leftforepaw_session_likelihood(neg_closestvalue) >= 0.95,'Error with DLC point processing')
        if neg_closestvalue-1 > 0
            assert(leftforepaw_session_likelihood(neg_closestvalue-1) < 0.95,'Error with DLC point processing')
        end
    end
    
    % Replace the bad outlier frame x and y values with mean of closest frame values
    leftforepaw_session_movement(:,bad_frame) = (leftforepaw_session_movement(:,neg_closestvalue) + leftforepaw_session_movement(:,pos_closestvalue))/2;
    
    clear neg_diff pos_only_neg_diff neg_closestvalue pos_diff pos_only_npos_diff pos_closestvalue
end

% Second step: Replace outliers
% Outlier removal Processing to remove outliers
time = [1:length(rightforepaw_session_movement)];
rightforepaw_session_movement(1,:) = filloutliers(rightforepaw_session_movement(1,:),'previous','movmedian',5,'SamplePoints',time); % remove outliers points that lie outside 3x mean absolute deviation from the local median within the 5 sampling points sliding window
rightforepaw_session_movement(2,:) = filloutliers(rightforepaw_session_movement(2,:),'previous','movmedian',5,'SamplePoints',time); % remove outliers points that lie outside 3x mean absolute deviation from the local median within the 5 sampling points sliding window
leftforepaw_session_movement(1,:) = filloutliers(leftforepaw_session_movement(1,:),'previous','movmedian',5,'SamplePoints',time); % remove outliers points that lie outside 3x mean absolute deviation from the local median within the 5 sampling points sliding window
leftforepaw_session_movement(2,:) = filloutliers(leftforepaw_session_movement(2,:),'previous','movmedian',5,'SamplePoints',time); % remove outliers points that lie outside 3x mean absolute deviation from the local median within the 5 sampling points sliding window

for trial_num = 1:length(trial_begin_img)
    
    % First find the average img bin size within the trial.
    trial_bin_size = round(mean(img_frame_end(trial_begin_img(trial_num): trial_end_img(trial_num)) - img_frame_begin(trial_begin_img(trial_num): trial_end_img(trial_num))));
    % This video bin size will be used to average all movie frames within bin
    video_bin_size = round(trial_bin_size / video_frame_rate);
    
    fs_wavesurfer_fs_movie_ratio = (img_frame_begin(trial_begin_img(trial_num)) / movie_trial_begin(trial_num));
    fs_wavesurfer_fs_movie_ratio_trial(trial_num) = (img_frame_begin(trial_begin_img(trial_num)) / movie_trial_begin(trial_num));
    
    movie_frame_begin_img = round(img_frame_begin(trial_begin_img(trial_num) : trial_end_img(trial_num)) ./ fs_wavesurfer_fs_movie_ratio);
    movie_frame_end_img = round(img_frame_end(trial_begin_img(trial_num) : trial_end_img(trial_num)) ./ fs_wavesurfer_fs_movie_ratio);
    
    % Start binning the x,y coordinates by taking mean of bin window
    for frame_num = 1:length(movie_frame_begin_img)
        binned_rightforepaw_movement{trial_num}(:,frame_num) = mean(rightforepaw_session_movement(:,movie_frame_begin_img(frame_num): movie_frame_end_img(frame_num)),2);
        binned_leftforepaw_movement{trial_num}(:,frame_num) = mean(leftforepaw_session_movement(:,movie_frame_begin_img(frame_num): movie_frame_end_img(frame_num)),2);
    end
    
    clear movie_frame_begin_img clear movie_frame_end_img
end

clear fs_wavesurfer_fs_movie_ratio clear movie_frame_begin_img movie_frame_end_img trial_bin_size video_bin_size

% After binning according to trial number, get movement between each frame
for trial_num = 1:length(movie_trial_begin)
    binned_rightforepaw_movement{trial_num} = binned_rightforepaw_movement{trial_num}';
    binned_leftforepaw_movement{trial_num} = binned_leftforepaw_movement{trial_num}';
    for pixel = 1:length(binned_rightforepaw_movement{trial_num})-1
        rightforepaw_euclidean{trial_num}(1,pixel) = pdist(binned_rightforepaw_movement{trial_num}(pixel:pixel+1,:),'euclidean'); % Get euclidean distance between time point t and t+1
        leftforepaw_euclidean{trial_num}(1,pixel) = pdist(binned_leftforepaw_movement{trial_num}(pixel:pixel+1,:),'euclidean'); % Get euclidean distance between time point t and t+1
    end
    
    % Pre-processing: Smoothing by 3 points/sliding window
    processed_rightforepaw_euclidean{trial_num} = smooth(rightforepaw_euclidean{trial_num},3)';
    processed_leftforepaw_euclidean{trial_num} = smooth(leftforepaw_euclidean{trial_num},3)';
end

% % Concatenate across trials without ITIs.
right_forepaw_euclidean_session_img = zeros(length(img_frame_begin),1)';
left_forepaw_euclidean_session_img = zeros(length(img_frame_begin),1)';

for trial_num = 1:length(trial_begin)
    right_forepaw_euclidean_session_img(:,trial_begin_img(trial_num): (trial_begin_img(trial_num)+length(processed_rightforepaw_euclidean{trial_num})-1) ) = processed_rightforepaw_euclidean{trial_num};
    left_forepaw_euclidean_session_img(:,trial_begin_img(trial_num): (trial_begin_img(trial_num)+length(processed_leftforepaw_euclidean{trial_num})-1) ) = processed_leftforepaw_euclidean{trial_num};
end

right_forepaw_euclidean_session_img = right_forepaw_euclidean_session_img';
left_forepaw_euclidean_session_img = left_forepaw_euclidean_session_img';

clear rightforepaw_movement left_forepaw_movement combined_csv_table csv_table sorted_csv_file_list name_position_vector csv_file_list listing start_position index idx pixel i ...
    csv_num

%% Generalized linear model: Initialise task variables (Tactile stimulus onset left, Tactile stimulus onset right,lick onset left, lick onset right, reward onset, delay onset, rightforepaw , leftforepaw).
tact_stim_onset_left = zeros(length(img_frame_begin),1);
tact_stim_onset_left(stim_middle_img(SessionData.TrialTypes == 1)) = 1;
tact_stim_onset_right = zeros(length(img_frame_begin),1);
tact_stim_onset_right(stim_middle_img(SessionData.TrialTypes == 2)) = 1;

% Reward onset is defined by first left or right lick, given that the trial type is correct
% Edited - push reward onset to 0.5s post lick onset
reward_onset_img_trial = nan(1,length(trial_type));
reward_onset_img_trial(find(trial_type == 1)) = post_go_cue_lick_left_img(find(trial_type == 1)) + round(fs_image*0.5);
reward_onset_img_trial(find(trial_type == 3)) = post_go_cue_lick_right_img(find(trial_type == 3)) + round(fs_image*0.5);
reward_onset_img_trial = reward_onset_img_trial(~isnan(reward_onset_img_trial)); % Remove nan
reward_onset = zeros(length(img_frame_begin),1); % Initialise reward onset vector
reward_onset(reward_onset_img_trial) = 1;

% post_go_cue_lick_center_img = post_go_cue_lick_center_img(~isnan(post_go_cue_lick_center_img)); % Remove nans.
post_go_cue_lick_left_img = post_go_cue_lick_left_img(~isnan(post_go_cue_lick_left_img)); % Remove nans.
post_go_cue_lick_right_img = post_go_cue_lick_right_img(~isnan(post_go_cue_lick_right_img)); % Remove nans.

% Using first lick for left or right to create predictors
lick_left_onset = zeros(length(img_frame_begin),1);
lick_left_onset(post_go_cue_lick_left_img) = 1;
lick_right_onset = zeros(length(img_frame_begin),1);
lick_right_onset(post_go_cue_lick_right_img) = 1;

tact_stim_onset_left = zeros(length(img_frame_begin),1);
tact_stim_onset_left(trial_begin_img(SessionData.TrialTypes == 1)) = 1;

% Left delay
delay_onset_left = zeros(length(img_frame_begin),1);
delay_onset_left(delay_begin_img(SessionData.TrialTypes == 1)) = 1;
delay_offset_left = zeros(length(img_frame_begin),1);
delay_offset_left(delay_end_img(SessionData.TrialTypes == 1)) = 1;
delay_middle_left = zeros(length(img_frame_begin),1);
delay_middle_left(delay_middle_img(SessionData.TrialTypes == 1)) = 1;

% Right delay
delay_onset_right = zeros(length(img_frame_begin),1);
delay_onset_right(delay_begin_img(SessionData.TrialTypes == 2)) = 1;
delay_offset_right = zeros(length(img_frame_begin),1);
delay_offset_right(delay_end_img(SessionData.TrialTypes == 2)) = 1;
delay_middle_right = zeros(length(img_frame_begin),1);
delay_middle_right(delay_middle_img(SessionData.TrialTypes == 2)) = 1;

go_cue_onset = zeros(length(img_frame_begin),1);
go_cue_onset(go_cue_begin_img,1) = 1;

%% Convolution of task variables with basis functions. This part is based on Pillow lab's code.
% Stimulus onset - 8 evenly spaced raised cosine basis function.
% Raised cosine basis function.
clear nkbins nBases ttb dbcenter width bcenters bfun BBstm

% Edited - duration of basis function
resp_basis_duration = 1.5; % Reward & tact stim & left and right lick
stim_basis_duration = 1;
reward_basis_duration = 1;
delay_basis_duration = 1.5;
movement_basis_duration = 1;

% Response
% Fix the distance between 2 cosine from peak to peak
fixed_dbcenter = 2.5;

% Generating basis functions for resp;
basis_duration = resp_basis_duration;
nkbins = round(basis_duration/(1/fs_image)); % Bin number over 6 s. % Over 6s how many img bins are there
nkbins = round((nkbins)/2)*2; % Round to the nearest even number.

% For raised cosine, the spacing between the centers mu2 1st be 1/4 of the width of the cosine.
width = ceil(4*fixed_dbcenter); % Width of each bump.

% Check how many basis functions are needed to fit within the nkbins given the fixed width
center = width/4; % Same as fixed_dbcenter
bcenter = (width/2-2) : center : (nkbins); % Manually fix the start of the basis function as -2
resp_nBases = length(bcenter);
ttb = repmat((1:nkbins)',1,resp_nBases); % Time indices for basis.
bfun = @(x,period)((abs(x/period) < 0.5).*(cos(x*2*pi/period)*0.5 + 0.5));
BBstm = bfun(ttb - repmat(bcenter,nkbins,1),width);
resp_nkbins = nkbins; % To save specific basis function for specific cues/delay/movement
resp_BBstm = BBstm;

clear nBases basis_duration nkbins ttb dbcenter width bcenters bfun BBstm

% Stimulus
% Fix the distance between 2 cosine
fixed_dbcenter = 2.5; clear nBases

% Generating basis functions for resp;
basis_duration = stim_basis_duration;
nkbins = round(basis_duration/(1/fs_image)); % Bin number over 6 s. % Over 6s how many img bins are there
nkbins = round((nkbins)/2)*2; % Round to the nearest even number.
width = ceil(4*fixed_dbcenter); % Width of each bump.

% Check how many basis functions are needed to fit within the nkbins given the fixed width
% For raised cosine, the spacing between the centers mu2 1st be 1/4 of the width of the cosine.
center = width/4;
bcenter = (width/2-2) : center : (nkbins);
stim_nBases = length(bcenter);
ttb = repmat((1:nkbins)',1,stim_nBases); % Time indices for basis.
bfun = @(x,period)((abs(x/period) < 0.5).*(cos(x*2*pi/period)*0.5 + 0.5));
BBstm = bfun(ttb - repmat(bcenter,nkbins,1),width);
stim_nkbins = nkbins; % To save specific basis function for specific cues/delay/movement
stim_BBstm = BBstm;
clear nBases basis_duration nkbins ttb dbcenter width bcenters bfun BBstm

% Delay
% Fix the distance between 2 cosine
fixed_dbcenter = 2.5; clear nBases

% Generating basis functions for resp;
basis_duration = delay_basis_duration;
nkbins = round(basis_duration/(1/fs_image)); % Bin number over 6 s. % Over 6s how many img bins are there
nkbins = round((nkbins)/2)*2; % Round to the nearest even number.
width = ceil(4*fixed_dbcenter); % Width of each bump.

% Check how many basis functions are needed to fit within the nkbins given the fixed width
% For raised cosine, the spacing between the centers mu2 1st be 1/4 of the width of the cosine.
center = width/4;
bcenter = (width/2-2) : center : (nkbins);
delay_nBases = length(bcenter);
ttb = repmat((1:nkbins)',1,delay_nBases); % Time indices for basis.
bfun = @(x,period)((abs(x/period) < 0.5).*(cos(x*2*pi/period)*0.5 + 0.5));
BBstm = bfun(ttb - repmat(bcenter,nkbins,1),width);
delay_nkbins = nkbins; % To save specific basis function for specific cues/delay/movement
delay_BBstm = BBstm;
clear nBases basis_duration nkbins ttb dbcenter width bcenters bfun BBstm

% reward
% Fix the distance between 2 cosine
fixed_dbcenter = 2.5; clear nBases

% Generating basis functions for resp;
basis_duration = reward_basis_duration;
nkbins = round(basis_duration/(1/fs_image)); % Bin number over 6 s. % Over 6s how many img bins are there
nkbins = round((nkbins)/2)*2; % Round to the nearest even number.
width = ceil(4*fixed_dbcenter); % Width of each bump.

% Check how many basis functions are needed to fit within the nkbins given the fixed width
% For raised cosine, the spacing between the centers mu2 1st be 1/4 of the width of the cosine.
center = width/4;
bcenter = (width/2-2) : center : (nkbins);
reward_nBases = length(bcenter);
ttb = repmat((1:nkbins)',1,reward_nBases); % Time indices for basis.
bfun = @(x,period)((abs(x/period) < 0.5).*(cos(x*2*pi/period)*0.5 + 0.5));
BBstm = bfun(ttb - repmat(bcenter,nkbins,1),width);
reward_nkbins = nkbins; % To save specific basis function for specific cues/reward/movement
reward_BBstm = BBstm;
clear nBases basis_duration nkbins ttb dbcenter width bcenters bfun BBstm

% movement
% Fix the distance between 2 cosine
fixed_dbcenter = 2.5; clear nBases

% Generating basis functions for resp;
basis_duration = movement_basis_duration;
nkbins = round(basis_duration/(1/fs_image)); % Bin number over 6 s. % Over 6s how many img bins are there
nkbins = round((nkbins)/2)*2; % Round to the nearest even number.
width = ceil(4*fixed_dbcenter); % Width of each bump.

% Check how many basis functions are needed to fit within the nkbins given the fixed width
% For raised cosine, the spacing between the centers mu2 1st be 1/4 of the width of the cosine.
center = width/4;
bcenter = (width/2-2) : center : (nkbins);
movement_nBases = length(bcenter);
ttb = repmat((1:nkbins)',1,movement_nBases); % Time indices for basis.
bfun = @(x,period)((abs(x/period) < 0.5).*(cos(x*2*pi/period)*0.5 + 0.5));
BBstm = bfun(ttb - repmat(bcenter,nkbins,1),width);
movement_nkbins = nkbins; % To save specific basis function for specific cues/movement/movement
movement_BBstm = BBstm;
clear nBases basis_duration nkbins ttb dbcenter width bcenters bfun BBstm

% Convolution (aligned to the center).
clear tact_stim_onset_left_convolved tact_stim_onset_right_convolved reward_onset_convolved delay_onset_convolved delay_offset_convolved delay_middle_convolved go_cue_convolved rightforepaw_movement_convolved leftforepaw_movement_convolved...
    left_lick_convolved right_lick_convolved delay_left_onset_convolved delay_left_offset_convolved delay_left_middle_convolved delay_right_onset_convolved delay_right_offset_convolved delay_right_middle_convolved...
    lick_left_onset_convolved lick_right_onset_convolved

% Stim
for frame_num = (stim_nkbins/2):(length(img_frame_end) - (stim_nkbins/2))
    tact_stim_onset_left_convolved(frame_num,:) = tact_stim_onset_left((frame_num - stim_nkbins/2 + 1):(frame_num + (stim_nkbins/2)))'*stim_BBstm;
    tact_stim_onset_right_convolved(frame_num,:) = tact_stim_onset_right((frame_num - stim_nkbins/2 + 1):(frame_num + (stim_nkbins/2)))'*stim_BBstm;
end

% Reward
for frame_num = (reward_nkbins/2):(length(img_frame_end) - (reward_nkbins/2))
    reward_onset_convolved(frame_num,:) = reward_onset((frame_num - reward_nkbins/2 + 1):(frame_num + (reward_nkbins/2)))'*reward_BBstm;
end

% Response
for frame_num = (resp_nkbins/2):(length(img_frame_end) - (resp_nkbins/2))
    lick_left_onset_convolved(frame_num,:) = lick_left_onset((frame_num - resp_nkbins/2 + 1):(frame_num + (resp_nkbins/2)))'*resp_BBstm;
    lick_right_onset_convolved(frame_num,:) = lick_right_onset((frame_num - resp_nkbins/2 + 1):(frame_num + (resp_nkbins/2)))'*resp_BBstm;
end

% Delay
for frame_num = (delay_nkbins/2):(length(img_frame_end) - (delay_nkbins/2))
    
    % Edited 5/5/21 to add trial specific delay
    delay_left_onset_convolved(frame_num,:) = delay_onset_left((frame_num - delay_nkbins/2 + 1):(frame_num + (delay_nkbins/2)))'*delay_BBstm;
    delay_left_offset_convolved(frame_num,:) = delay_offset_left((frame_num - delay_nkbins/2 + 1):(frame_num + (delay_nkbins/2)))'*delay_BBstm;
    delay_left_middle_convolved(frame_num,:) = delay_middle_left((frame_num - delay_nkbins/2 + 1):(frame_num + (delay_nkbins/2)))'*delay_BBstm;
    
    delay_right_onset_convolved(frame_num,:) = delay_onset_right((frame_num - delay_nkbins/2 + 1):(frame_num + (delay_nkbins/2)))'*delay_BBstm;
    delay_right_offset_convolved(frame_num,:) = delay_offset_right((frame_num - delay_nkbins/2 + 1):(frame_num + (delay_nkbins/2)))'*delay_BBstm;
    delay_right_middle_convolved(frame_num,:) = delay_middle_right((frame_num - delay_nkbins/2 + 1):(frame_num + (delay_nkbins/2)))'*delay_BBstm;
end

% Movement
for frame_num = (movement_nkbins/2):(length(img_frame_end) - (movement_nkbins/2))
    rightforepaw_movement_convolved(frame_num,:) = right_forepaw_euclidean_session_img((frame_num - movement_nkbins/2 + 1):(frame_num + (movement_nkbins/2)))'*movement_BBstm;
    leftforepaw_movement_convolved(frame_num,:) = left_forepaw_euclidean_session_img((frame_num - movement_nkbins/2 + 1):(frame_num + (movement_nkbins/2)))'*movement_BBstm;
end

%% Normalize each variable by zscore.
% Stimulus onset.
tact_stim_onset_left_convolved_all = tact_stim_onset_left_convolved(:);% Combining all 8 basis functions to obtain z score
norm_tact_stim_onset_left_convolved_all = zscore(tact_stim_onset_left_convolved_all); % Transforming to zscore in order to compare across all task variables
norm_tact_stim_onset_left_convolved = reshape(norm_tact_stim_onset_left_convolved_all,[size(tact_stim_onset_left_convolved,1),size(tact_stim_onset_left_convolved,2)]); % Reshaping it back to the 8 basis function
tact_stim_onset_right_convolved_all = tact_stim_onset_right_convolved(:);
norm_tact_stim_onset_right_convolved_all = zscore(tact_stim_onset_right_convolved_all);
norm_tact_stim_onset_right_convolved = reshape(norm_tact_stim_onset_right_convolved_all,[size(tact_stim_onset_right_convolved,1),size(tact_stim_onset_right_convolved,2)]);

% Reward onset.
reward_onset_convolved_all = reward_onset_convolved(:);
norm_reward_onset_convolved_all = zscore(reward_onset_convolved_all);
norm_reward_onset_convolved = reshape(norm_reward_onset_convolved_all,[size(reward_onset_convolved,1),size(reward_onset_convolved,2)]);
assert(isequal(norm_reward_onset_convolved(500:1000,1),norm_reward_onset_convolved_all(500:1000)),'reshaping is wrong')

% Left lick onset
lick_left_onset_convolved_all = lick_left_onset_convolved(:);
norm_lick_left_onset_convolved_all = zscore(lick_left_onset_convolved_all);
norm_lick_left_convolved = reshape(norm_lick_left_onset_convolved_all,[size(lick_left_onset_convolved,1),size(lick_left_onset_convolved,2)]);

% Right lick onset
lick_right_onset_convolved_all = lick_right_onset_convolved(:);
norm_lick_right_onset_convolved_all = zscore(lick_right_onset_convolved_all);
norm_lick_right_convolved = reshape(norm_lick_right_onset_convolved_all,[size(lick_right_onset_convolved,1),size(lick_right_onset_convolved,2)]);

% Delay onset, offset and middle
delay_left_onset_convolved_all = delay_left_onset_convolved(:);
norm_delay_left_onset_convolved_all = zscore(delay_left_onset_convolved_all);
norm_delay_left_onset_convolved = reshape(norm_delay_left_onset_convolved_all,[size(delay_left_onset_convolved,1),size(delay_left_onset_convolved,2)]);

delay_left_offset_convolved_all = delay_left_offset_convolved(:);
norm_delay_left_offset_convolved_all = zscore(delay_left_offset_convolved_all);
norm_delay_left_offset_convolved = reshape(norm_delay_left_offset_convolved_all,[size(delay_left_offset_convolved,1),size(delay_left_offset_convolved,2)]);

delay_left_middle_convolved_all = delay_left_middle_convolved(:);
norm_delay_left_middle_convolved_all = zscore(delay_left_middle_convolved_all);
norm_delay_left_middle_convolved = reshape(norm_delay_left_middle_convolved_all,[size(delay_left_middle_convolved,1),size(delay_left_middle_convolved,2)]);

% Delay onset, offset and middle
delay_right_onset_convolved_all = delay_right_onset_convolved(:);
norm_delay_right_onset_convolved_all = zscore(delay_right_onset_convolved_all);
norm_delay_right_onset_convolved = reshape(norm_delay_right_onset_convolved_all,[size(delay_right_onset_convolved,1),size(delay_right_onset_convolved,2)]);

delay_right_offset_convolved_all = delay_right_offset_convolved(:);
norm_delay_right_offset_convolved_all = zscore(delay_right_offset_convolved_all);
norm_delay_right_offset_convolved = reshape(norm_delay_right_offset_convolved_all,[size(delay_right_offset_convolved,1),size(delay_right_offset_convolved,2)]);

delay_right_middle_convolved_all = delay_right_middle_convolved(:);
norm_delay_right_middle_convolved_all = zscore(delay_right_middle_convolved_all);
norm_delay_right_middle_convolved = reshape(norm_delay_right_middle_convolved_all,[size(delay_right_middle_convolved,1),size(delay_right_middle_convolved,2)]);

% Forepaw movement
rightforepaw_movement_convolved_all = rightforepaw_movement_convolved(:);
norm_rightforepaw_movement_convolved_all = zscore(rightforepaw_movement_convolved_all);
norm_rightforepaw_movement_convolved = reshape(norm_rightforepaw_movement_convolved_all,[size(rightforepaw_movement_convolved,1),size(rightforepaw_movement_convolved,2)]);

leftforepaw_movement_convolved_all = leftforepaw_movement_convolved(:);
norm_leftforepaw_movement_convolved_all = zscore(leftforepaw_movement_convolved_all);
norm_leftforepaw_movement_convolved = reshape(norm_leftforepaw_movement_convolved_all,[size(leftforepaw_movement_convolved,1),size(leftforepaw_movement_convolved,2)]);

% Edited 5/5/21 to add delay_left or delay_right
% Divide into different trial types and get trial sections (-4 second from trial onset to ITI onset).
for trial_num = 1:length(trial_begin)
    norm_tact_stim_onset_left_convolved_trial{trial_num} = norm_tact_stim_onset_left_convolved((trial_begin_img(trial_num) - 4.*round(fs_image)):trial_end_img(trial_num),:);
    norm_tact_stim_onset_right_convolved_trial{trial_num} = norm_tact_stim_onset_right_convolved((trial_begin_img(trial_num) - 4.*round(fs_image)):trial_end_img(trial_num),:);
    norm_delay_left_middle_convolved_trial{trial_num} = norm_delay_left_middle_convolved((trial_begin_img(trial_num) - 4.*round(fs_image)):trial_end_img(trial_num),:);
    norm_delay_right_middle_convolved_trial{trial_num} = norm_delay_right_middle_convolved((trial_begin_img(trial_num) - 4.*round(fs_image)):trial_end_img(trial_num),:);
    
    norm_lick_right_onset_convolved_trial{trial_num} = norm_lick_right_convolved((trial_begin_img(trial_num) - 4.*round(fs_image)):trial_end_img(trial_num),:);
    norm_lick_left_onset_convolved_trial{trial_num} = norm_lick_left_convolved((trial_begin_img(trial_num) - 4.*round(fs_image)):trial_end_img(trial_num),:);
    norm_reward_onset_convolved_trial{trial_num} = norm_reward_onset_convolved((trial_begin_img(trial_num) - 4.*round(fs_image)):trial_end_img(trial_num),:);
    norm_rightforepaw_movement_convolved_trial{trial_num} = norm_rightforepaw_movement_convolved((trial_begin_img(trial_num) - 4.*round(fs_image)):trial_end_img(trial_num),:);
    norm_leftforepaw_movement_convolved_trial{trial_num} = norm_leftforepaw_movement_convolved((trial_begin_img(trial_num) - 4.*round(fs_image)):trial_end_img(trial_num),:);
end

%% Sort spike data into trial types and training/test data
% Do not consider non-active neurons.
active_cell_threshold = 0; % Adjust threshold accordingly to obtain active cells
for region_num = 1:region
    resp_var{region_num} = [];
    if size(spike_all{region_num},1) >= 1 % Only consider regions with more than 1 cell
        for cell_num = 1:size(spike_all{region_num},1)
            smoothed_spike_all{region_num}(cell_num,:) = smooth(spike_all{region_num}(cell_num,1:trial_end_img(end)+2)',round(fs_image./2)); % Smooth and add 2 frames to the end for coupling purposes
        end
        active_cell{region_num} = find(max(smoothed_spike_all{region_num},[],2) >= active_cell_threshold); % Only consider cells with max activity >= 5.
        resp_var{region_num} = smoothed_spike_all{region_num}(active_cell{region_num},:);
        
    end
end

% Note that additional resp_var_trial_coupling was added for frame shift reasons
for region_num = 1:region
    if size(resp_var{region_num},1) >= 1 % Only consider regions with more than 1 cell
        for trial_num = 1:length(trial_begin)
            resp_var_trial{region_num}{trial_num} = resp_var{region_num}(:,(trial_begin_img(trial_num) - 4.*round(fs_image)):trial_end_img(trial_num)); % Takes into account the 4s of the previous ITI
            resp_var_trial_coupling{region_num}{trial_num} = resp_var{region_num}(:,(trial_begin_img(trial_num)-2 - (4.*round(fs_image)):trial_end_img(trial_num))); % 2 additional frames added to the end to obtain t1 and t2 shifts
        end
    end
end

% Concatenate across regions.
for trial_num = 1:length(trial_begin)
    resp_var_trial_concat{trial_num} = [];
    resp_var_trial_concat_coupling{trial_num} = [];
    
    for region_num = 1:region
        if size(resp_var{region_num},1) >= 1
            resp_var_trial_concat{trial_num} = [resp_var_trial_concat{trial_num},resp_var_trial{region_num}{trial_num}'];
            resp_var_trial_concat_coupling{trial_num} = [resp_var_trial_concat_coupling{trial_num},resp_var_trial_coupling{region_num}{trial_num}'];
        end
    end
end

% Stratifying trial types (correct left, incorrect left, etc) into training or test dataset
seed = 1;
rng(seed); % To be deterministic.

trial_type_1_idx = find(trial_type == 1); % Correct left
trial_type_2_idx = find(trial_type == 2); % Incorect Left
trial_type_3_idx = find(trial_type == 3); % Correct Right
trial_type_4_idx = find(trial_type == 4); % Incorrect Right
trial_type_0_idx = find(trial_type == 0); % No response

trial_type_1_train = randsample(trial_type_1_idx,floor(length(trial_type_1_idx)*0.7),false);
trial_type_1_test = trial_type_1_idx(~ismember(trial_type_1_idx,trial_type_1_train));
trial_type_2_train = randsample(trial_type_2_idx,floor(length(trial_type_2_idx)*0.7),false);
trial_type_2_test = trial_type_2_idx(~ismember(trial_type_2_idx,trial_type_2_train));
trial_type_3_train = randsample(trial_type_3_idx,floor(length(trial_type_3_idx)*0.7),false);
trial_type_3_test = trial_type_3_idx(~ismember(trial_type_3_idx,trial_type_3_train));
trial_type_4_train = randsample(trial_type_4_idx,floor(length(trial_type_4_idx)*0.7),false);
trial_type_4_test = trial_type_4_idx(~ismember(trial_type_4_idx,trial_type_4_train));
trial_type_0_train = randsample(trial_type_0_idx,floor(length(trial_type_0_idx)*0.7),false);
trial_type_0_test = trial_type_0_idx(~ismember(trial_type_0_idx,trial_type_0_train));

% Check if assignment of trial types encompasses all 180 trials
assert(isequal(sort([trial_type_2_train;trial_type_2_test]),trial_type_2_idx),'Trial type error')
assert(isequal(sort([trial_type_1_train;trial_type_1_test;trial_type_2_train;trial_type_2_test;trial_type_3_train;trial_type_3_test;trial_type_4_train;...
    trial_type_4_test;trial_type_0_train;trial_type_0_test]),[1:180]'),'Trial type error')

task1_train_idx = [trial_type_1_train; trial_type_2_train; trial_type_3_train; trial_type_4_train; trial_type_0_train]; % Concatenate all the training data from trial type 1 and 2
task1_train_idx = sort(task1_train_idx);
task1_test_idx = [trial_type_1_test; trial_type_2_test; trial_type_3_test; trial_type_4_test; trial_type_0_test];
task1_test_idx = sort(task1_test_idx);
assert(isequal(sort([task1_train_idx;task1_test_idx]),[1:180]'),'Trial type error')

% Edited to remove go-cue
tact_stim_onset_left_idx = 1:stim_nBases;
tact_stim_onset_right_idx = (max(tact_stim_onset_left_idx)+1) : (max(tact_stim_onset_left_idx)+stim_nBases);
delay_left_middle_idx = (max(tact_stim_onset_right_idx)+1) : (max(tact_stim_onset_right_idx)+delay_nBases);
delay_right_middle_idx = (max(delay_left_middle_idx)+1) : (max(delay_left_middle_idx)+delay_nBases);
% go_cue_idx = 17:18;
left_lick_onset_idx = (max(delay_right_middle_idx)+1) : (max(delay_right_middle_idx)+resp_nBases);
right_lick_onset_idx = (max(left_lick_onset_idx)+1) : (max(left_lick_onset_idx)+resp_nBases);
reward_onset_idx = (max(right_lick_onset_idx)+1) : (max(right_lick_onset_idx)+reward_nBases);
leftforepaw_movement_idx = (max(reward_onset_idx)+1) : (max(reward_onset_idx)+movement_nBases);
rightforepaw_movement_idx = (max(leftforepaw_movement_idx)+1) : (max(leftforepaw_movement_idx)+movement_nBases);

% Training dataset.
clear task1_var_train_temp
for trial_num = 1:length(task1_train_idx)
    task1_var_train_temp{trial_num} = [
        norm_tact_stim_onset_left_convolved_trial{task1_train_idx(trial_num)},...
        norm_tact_stim_onset_right_convolved_trial{task1_train_idx(trial_num)},...
        norm_delay_left_middle_convolved_trial{task1_train_idx(trial_num)},...
        norm_delay_right_middle_convolved_trial{task1_train_idx(trial_num)},...
        norm_lick_left_onset_convolved_trial{task1_train_idx(trial_num)},...
        norm_lick_right_onset_convolved_trial{task1_train_idx(trial_num)},...
        norm_reward_onset_convolved_trial{task1_train_idx(trial_num)},...
        norm_leftforepaw_movement_convolved_trial{task1_train_idx(trial_num)},...
        norm_rightforepaw_movement_convolved_trial{task1_train_idx(trial_num)}];
end

assert(size(task1_var_train_temp{1},2)==max(rightforepaw_movement_idx))

% Concatenate all training variables across the entire session.
task1_var_train = [];
resp_task1_var_train = [];
for trial_num = 1:length(task1_train_idx)
    task1_var_train = [task1_var_train;task1_var_train_temp{trial_num}];
    resp_task1_var_train = [resp_task1_var_train;resp_var_trial_concat{task1_train_idx(trial_num)}];
end

clear task1_var_test_temp
% Test dataset. (Repeating the same for what we did for training dataset)
for trial_num = 1:length(task1_test_idx)
    task1_var_test_temp{trial_num} = [
        norm_tact_stim_onset_left_convolved_trial{task1_test_idx(trial_num)},...
        norm_tact_stim_onset_right_convolved_trial{task1_test_idx(trial_num)},...
        norm_delay_left_middle_convolved_trial{task1_test_idx(trial_num)},...
        norm_delay_right_middle_convolved_trial{task1_test_idx(trial_num)},...
        norm_lick_left_onset_convolved_trial{task1_test_idx(trial_num)},...
        norm_lick_right_onset_convolved_trial{task1_test_idx(trial_num)},...
        norm_reward_onset_convolved_trial{task1_test_idx(trial_num)},...
        norm_leftforepaw_movement_convolved_trial{task1_test_idx(trial_num)},...
        norm_rightforepaw_movement_convolved_trial{task1_test_idx(trial_num)}];
end
assert(size(task1_var_test_temp{1},2)==max(rightforepaw_movement_idx))

% Concatenate.
task1_var_test = [];
resp_task1_var_test = [];
for trial_num = 1:length(task1_test_idx)
    task1_var_test = [task1_var_test;task1_var_test_temp{trial_num}];
    resp_task1_var_test = [resp_task1_var_test;resp_var_trial_concat{task1_test_idx(trial_num)}];
end

%% Coupling between neurons (Normalise variables -> sort into test/train)
% Create new handle for task variables with coupling predictors
task1_var_coupling_train = task1_var_train;
task1_var_coupling_test = task1_var_test;

for trial_num = 1 : size(resp_var_trial_concat_coupling,2)
    % Edited 14 Jul 21 to change the direction of predictor
    resp_var_trial_concat_coupling_t1{trial_num} = resp_var_trial_concat_coupling{trial_num}(2:end-1,:); % The additional 2 frames that was added is removed in this process
    resp_var_trial_concat_coupling_t2{trial_num} = resp_var_trial_concat_coupling{trial_num}(1:end-2,:);
    
end

% Concatenate across all trials for normalising
resp_task1_var_t1 = [];
resp_task1_var_t2 = [];
for trial_num = 1: size(resp_var_trial_concat_coupling,2)
    resp_task1_var_t1 = [resp_task1_var_t1 ; resp_var_trial_concat_coupling_t1{trial_num}];
    resp_task1_var_t2 = [resp_task1_var_t2 ; resp_var_trial_concat_coupling_t2{trial_num}];
end

% Normalise predictors across all cells
resp_task1_var_t1_cat = resp_task1_var_t1(:);
norm_resp_task1_var_t1_cat = zscore(resp_task1_var_t1_cat);
norm_resp_task1_var_t1 = reshape(norm_resp_task1_var_t1_cat,[size(resp_task1_var_t1,1),size(resp_task1_var_t1,2)]);
assert(isequal(norm_resp_task1_var_t1_cat(1:100),norm_resp_task1_var_t1(1:100,1)),'reshaping is wrong')

resp_task1_var_t2_cat = resp_task1_var_t2(:);
norm_resp_task1_var_t2_cat = zscore(resp_task1_var_t2_cat);
norm_resp_task1_var_t2 = reshape(norm_resp_task1_var_t2_cat,[size(resp_task1_var_t2,1),size(resp_task1_var_t2,2)]);

% Rearrange back into trials after normalising to sort into training and test dataset
start = 0;
for trial_num = 1 : size(resp_var_trial_concat_coupling,2)
    norm_resp_task1_var_t1_trial{trial_num} = norm_resp_task1_var_t1( start+1 : (start+size(resp_var_trial_concat_coupling_t1{trial_num},1)),:);
    norm_resp_task1_var_t2_trial{trial_num} = norm_resp_task1_var_t2( start+1 : (start+size(resp_var_trial_concat_coupling_t2{trial_num},1)),:);
    
    start = start + size(resp_var_trial_concat_coupling_t1{trial_num},1);
end

% Sort training set
norm_resp_task1_var_t1_train = [];
norm_resp_task1_var_t2_train = [];
for i = 1:size(task1_train_idx,1)
    train_idx = task1_train_idx(i);
    norm_resp_task1_var_t1_train = [norm_resp_task1_var_t1_train;norm_resp_task1_var_t1_trial{train_idx}];
    norm_resp_task1_var_t2_train = [norm_resp_task1_var_t2_train;norm_resp_task1_var_t2_trial{train_idx}];
end
% Sort test set
norm_resp_task1_var_t1_test = [];
norm_resp_task1_var_t2_test = [];
for i = 1:size(task1_test_idx,1)
    test_idx = task1_test_idx(i);
    norm_resp_task1_var_t1_test = [norm_resp_task1_var_t1_test;norm_resp_task1_var_t1_trial{test_idx}];
    norm_resp_task1_var_t2_test = [norm_resp_task1_var_t2_test;norm_resp_task1_var_t2_trial{test_idx}];
end

% Arrange predictors t1 and t2 for the same cell side by side
norm_resp_task1_var_t1_t2_train = [];
predictor_labels = nan(size(task1_var_train,2),1); % Create predictor labels to keep track of stimulus and cell predictors
for cell_num = 1:size(norm_resp_task1_var_t1_train,2)
    norm_resp_task1_var_t1_t2_train = [norm_resp_task1_var_t1_t2_train norm_resp_task1_var_t1_train(:,cell_num) norm_resp_task1_var_t2_train(:,cell_num)];
    predictor_labels = [predictor_labels; cell_num; cell_num]; % Add cell number twice for t1 and t2
end

norm_resp_task1_var_t1_t2_test = [];
for cell_num = 1:size(norm_resp_task1_var_t1_test,2)
    norm_resp_task1_var_t1_t2_test = [norm_resp_task1_var_t1_t2_test norm_resp_task1_var_t1_test(:,cell_num) norm_resp_task1_var_t2_test(:,cell_num)];
end

% Add coupling predictors to the behavioural predictors
task1_var_coupling_train = [task1_var_coupling_train  norm_resp_task1_var_t1_t2_train];
task1_var_coupling_test = [task1_var_coupling_test  norm_resp_task1_var_t1_t2_test];

assert(size(task1_var_coupling_train,2) == (size(region_labels,1)*2 + size(task1_var_train,2)),'Number of predictors is wrong')
[a,~] = cellfun(@size,resp_var); % To test size
assert(size(task1_var_coupling_test,2) == (sum(a)*2 + size(task1_var_test_temp{1},2)),'Number of predictors is wrong')
clear a

disp('Saving')
cd(parameters.root_folder);

% Save all variables.
mkdir(parameters.root_folder,'results')
cd([parameters.root_folder,'/results'])
save(['GLM_parameters.mat'],'-v7.3')