%% Analyze Behaviour for DREADD animals
function [CNO_correct_rate, Saline_correct_rate] = analyze_behaviour_dreadd(DREADD_data)
for animal_num = 1:size(DREADD_data.cno.SessionData,1)
    for file_num = 1:size(DREADD_data.cno.SessionData,2)
        % Load Session data for both CNO and Saline
        CNO_SessionData = DREADD_data.cno.SessionData(animal_num,file_num);
        clear SessionData
        
        Saline_SessionData = DREADD_data.saline.SessionData(animal_num,file_num);
        clear SessionData
        
        % Get correct rate
        CNO_correct_rate(file_num,1,animal_num) = sum(CNO_SessionData.Outcomes(find(CNO_SessionData.TrialTypes <= 2))) / length(CNO_SessionData.Outcomes(find(CNO_SessionData.TrialTypes <= 2)));
        CNO_correct_rate(file_num,2,animal_num) = sum(CNO_SessionData.Outcomes(find(CNO_SessionData.TrialTypes >= 3))) / length(CNO_SessionData.Outcomes(find(CNO_SessionData.TrialTypes >= 3)));
        
        % Get correct rate
        Saline_correct_rate(file_num,1,animal_num) = sum(Saline_SessionData.Outcomes(find(Saline_SessionData.TrialTypes <= 2))) / length(Saline_SessionData.Outcomes(find(Saline_SessionData.TrialTypes <= 2)));
        Saline_correct_rate(file_num,2,animal_num) = sum(Saline_SessionData.Outcomes(find(Saline_SessionData.TrialTypes >= 3))) / length(Saline_SessionData.Outcomes(find(Saline_SessionData.TrialTypes >= 3)));
    end
end
end