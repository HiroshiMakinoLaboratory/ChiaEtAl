%% Run functions to analyse example trained RNN
% Select folder where 'RNN_data' folder is saved 
folder_path = uigetdir;

% Run demo function to get choice activity form demo data
example_RNN = get_example_choice_mode_RNN(folder_path, 'Distractor_trials');

% Run to plot population activity projected on choice axis for +PPC and +vS1
plot_RNN_choice_mode(example_RNN,1)