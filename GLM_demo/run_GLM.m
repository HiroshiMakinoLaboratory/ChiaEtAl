%% (Demo) Randomly separate the data into training dataset and test dataset and perform GLM.
% Load.
parameters.root_folder = uigetdir;
cd([parameters.root_folder,'/GLM_demo'])
load 'GLM_parameters.mat'

% Image parameters
fs_image = 9.352;
x_size = 500;
y_size = 500;
region = 8;
xy_pixel_ratio = 0.4;
alpha = 0.95; % For elastic net regularisation

% Assign indices for cvpartition manually.
n = [1:size(resp_task1_var_train,1)];
k_fold = 5;
cvp = cvpartition(n,'KFold',k_fold);
for k_fold_num = 1:k_fold
    cvp_idx{k_fold_num} = (round(size(resp_task1_var_train,1)./k_fold)*(k_fold_num - 1) + 1):round(size(resp_task1_var_train,1)./k_fold)*(k_fold_num);
    cvp.Impl.indices(cvp_idx{k_fold_num}) = k_fold_num;
end
cvp.Impl.indices = cvp.Impl.indices(1:length(n));

% Elastic net regularization for GLM. Alpha is a weight of lasso (L1) versus ridge (L2) optimization. Alpha = 1 mean lasso regulation.
for cell_num = 1 : size(resp_task1_var_train,2)
    disp(['Cell number: ',num2str(cell_num),' out of ', num2str(size(resp_task1_var_train,2)),])
    clear B FitInfo idxLambdaMinDeviance
    
    % Remove same cell to prevent autocorrelation
    same_cell_predictor_idx = find(predictor_labels == cell_num);
    
    % Replace the predictor time series by its mean
    task1_var_coupling_train_cells{cell_num} = task1_var_coupling_train;
    task1_var_coupling_train_cells{cell_num}(:,same_cell_predictor_idx) = repmat(mean(task1_var_coupling_train_cells{cell_num}(:,same_cell_predictor_idx)),size(task1_var_coupling_train_cells{cell_num}(:,same_cell_predictor_idx),1),1);
    
    % Replace the predictor time series by its mean
    task1_var_coupling_test_cells{cell_num} = task1_var_coupling_test;
    task1_var_coupling_test_cells{cell_num}(:,same_cell_predictor_idx) = repmat(mean(task1_var_coupling_test_cells{cell_num}(:,same_cell_predictor_idx)),size(task1_var_coupling_test_cells{cell_num}(:,same_cell_predictor_idx),1),1);
    
    clear same_cell_predictor_idx
    
    % B is the max likelihood fitted coefficient where row : predictors for 100 (col) values given different regularisation coefficient
    [B,FitInfo] = lassoglm(task1_var_coupling_train_cells{cell_num},resp_task1_var_train(:,cell_num),'poisson','Alpha',alpha,'CV',cvp,'Link','log','MaxIter',1e4,'NumLambda',100,'Options',statset('UseParallel',true,'Streams',RandStream('mlfg6331_64'),'UseSubstreams',true),'Standardize',false);
    
    % Lambda min deviance is the regularisation coefficient that gives the lowest deviance fit
    idxLambdaMinDeviance = FitInfo.IndexMinDeviance;
    B0(cell_num) = FitInfo.Intercept(idxLambdaMinDeviance); % For the lambda with min deviance, what is the fit intercept
    % Selecting the task predictors with the min deviance
    coeff(cell_num,:) = B(:,idxLambdaMinDeviance); % For the predictors with the min deviance, what is the coefficient for all predictors (colm
    predict_train(cell_num,:) = exp(B0(cell_num) + task1_var_coupling_train_cells{cell_num}*coeff(cell_num,:)'); % Get dot product of training task predictors and task coefficients
    predict_test(cell_num,:) = exp(B0(cell_num) + task1_var_coupling_test_cells{cell_num}*coeff(cell_num,:)');
    
    % Model assessment by measuing explained deviance based on Benjamin et al 2018.
    y_test(cell_num,:) = resp_task1_var_test(:,cell_num); % Actual neuronal activity
    y_hat_test(cell_num,:) = predict_test(cell_num,:); % Predicted neuronal activity - predicted using the task coefficients for each task predictors
    y_null_test(cell_num) = mean(resp_task1_var_test(:,cell_num));% Using mean of response as the null distribution)
    L1_test(cell_num) = sum(y_test(cell_num,:).*log(eps + y_hat_test(cell_num,:)) - y_hat_test(cell_num,:)); % Predicted model.
    L0_test(cell_num) = sum(y_test(cell_num,:).*log(eps + y_null_test(cell_num)) - y_null_test(cell_num)); % Null model.
    LS_test(cell_num) = sum(y_test(cell_num,:).*log(eps + y_test(cell_num,:)) - y_test(cell_num,:)); % Saturated model.
    explained_deviance_test(cell_num) = 1 - (LS_test(cell_num) - L1_test(cell_num))/(LS_test(cell_num) - L0_test(cell_num)); % Explained deviance using all task predictors .* all task coefficients
    
    % Put the result in the structure.
    GLM_result.cvp.Impl.indices = cvp.Impl.indices;
    GLM_result.B0(cell_num) = B0(cell_num);
    GLM_result.coeff(cell_num,:) = coeff(cell_num,:);
    GLM_result.predict_train(cell_num,:) = predict_train(cell_num,:);
    GLM_result.predict_test(cell_num,:) = predict_test(cell_num,:);
    GLM_result.y_test(cell_num,:) = y_test(cell_num,:);
    GLM_result.y_hat_test(cell_num,:) = y_hat_test(cell_num,:);
    GLM_result.y_null_test(cell_num) = y_null_test(cell_num);
    GLM_result.L1_test(cell_num) = L1_test(cell_num);
    GLM_result.L0_test(cell_num) = L0_test(cell_num);
    GLM_result.LS_test(cell_num) = LS_test(cell_num);
    GLM_result.explained_deviance_test(cell_num) = explained_deviance_test(cell_num);
    
    % Save GLM results.
    save(['GLM_result.mat'],'GLM_result')
end