%% Analyze Behaviour for Optogenetics animals
function [ppc_correct_rate, vs1_correct_rate] = analyze_behaviour_opto(opto_data)
for animal_num = 1:size(opto_data.ppc.SessionData,1)
    for file_num = 1:size(opto_data.ppc.SessionData,2)
        
        % Load Session data for both CNO and Saline
        ppc_SessionData = opto_data.ppc.SessionData(animal_num,file_num);
        vs1_SessionData = opto_data.vs1.SessionData(animal_num,file_num);
        
        %% Get correct rate PPC
        % Control - no opto
        trial_12_idx =  find(ppc_SessionData.TrialTypes == 1 | ppc_SessionData.TrialTypes == 2);
        trial_12_idx(find(ppc_SessionData.Outcomes(trial_12_idx) >1)) = [];
        ppc_correct_rate(file_num,1,animal_num) = sum(ppc_SessionData.Outcomes(trial_12_idx)) / length(trial_12_idx);
        % Control - opto
        trial_34_idx =  find(ppc_SessionData.TrialTypes == 3 | ppc_SessionData.TrialTypes == 4);
        trial_34_idx(find(ppc_SessionData.Outcomes(trial_34_idx) >1)) = [];
        ppc_correct_rate(file_num,2,animal_num) = sum(ppc_SessionData.Outcomes(trial_34_idx)) / length(trial_34_idx);
        % Distractor - no opto
        trial_56_idx =  find(ppc_SessionData.TrialTypes == 5 | ppc_SessionData.TrialTypes == 6);
        trial_56_idx(find(ppc_SessionData.Outcomes(trial_56_idx) >1)) = [];
        ppc_correct_rate(file_num,3,animal_num) = sum(ppc_SessionData.Outcomes(trial_56_idx)) / length(trial_56_idx);
        % Distractor - opto
        trial_78_idx =  find(ppc_SessionData.TrialTypes == 7 | ppc_SessionData.TrialTypes == 8);
        trial_78_idx(find(ppc_SessionData.Outcomes(trial_78_idx) >1)) = [];
        ppc_correct_rate(file_num,4,animal_num) = sum(ppc_SessionData.Outcomes(trial_78_idx)) / length(trial_78_idx);
        
        %% Get correct rate VS1
        trial_12_idx =  find(vs1_SessionData.TrialTypes == 1 | vs1_SessionData.TrialTypes == 2);
        trial_12_idx(find(vs1_SessionData.Outcomes(trial_12_idx) >1)) = [];
        vs1_correct_rate(file_num,1,animal_num) = sum(vs1_SessionData.Outcomes(trial_12_idx)) / length(trial_12_idx);
        
        trial_34_idx =  find(vs1_SessionData.TrialTypes == 3 | vs1_SessionData.TrialTypes == 4);
        trial_34_idx(find(vs1_SessionData.Outcomes(trial_34_idx) >1)) = [];
        vs1_correct_rate(file_num,2,animal_num) = sum(vs1_SessionData.Outcomes(trial_34_idx)) / length(trial_34_idx);
        
        trial_56_idx =  find(vs1_SessionData.TrialTypes == 5 | vs1_SessionData.TrialTypes == 6);
        trial_56_idx(find(vs1_SessionData.Outcomes(trial_56_idx) >1)) = [];
        vs1_correct_rate(file_num,3,animal_num) = sum(vs1_SessionData.Outcomes(trial_56_idx)) / length(trial_56_idx);
        
        trial_78_idx =  find(vs1_SessionData.TrialTypes == 7 | vs1_SessionData.TrialTypes == 8);
        trial_78_idx(find(vs1_SessionData.Outcomes(trial_78_idx) >1)) = [];
        vs1_correct_rate(file_num,4,animal_num) = sum(vs1_SessionData.Outcomes(trial_78_idx)) / length(trial_78_idx);
    end
end
end