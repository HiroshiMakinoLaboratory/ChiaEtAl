%% Plot correlation between global enrichment indices and AUC
function plot_AUC_correlation(data,enrichment)
% Initialize data
enrich_idx = enrichment.enrich_idx;
sessions_color_map = [222 222 222; 145 145 145 ; 52 52 52]./255; % Color map

% Concat left and right enrichment indices together.
for session_num = 1:3
    sensory_enrich_idx{session_num} = [enrich_idx{1}{session_num}(~isnan(enrich_idx{1}{session_num})) enrich_idx{3}{session_num}(~isnan(enrich_idx{3}{session_num}))];
    delay_enrich_idx{session_num} = [enrich_idx{2}{session_num}(~isnan(enrich_idx{2}{session_num})) enrich_idx{4}{session_num}(~isnan(enrich_idx{4}{session_num}))];
end

%% Plot stim
figure('Position',[200,200,400,170],'Color','white','DefaultAxesFontSize',14);
y_lim = [0 0.8];x_lim = [0 1];
subplot(1,2,1)
for session_num = [1 3]
    
    AUC_data = data.AUC.choice_select_concat_stim_left_trial_auc{session_num};
    y = AUC_data(~isnan(AUC_data)); % Remove nans
    
    % Ensure that the length of AUC data and enrichment indices are the same
    assert(length(AUC_data) == length(sensory_enrich_idx{session_num}(~isnan(sensory_enrich_idx{session_num}))));
    x = sensory_enrich_idx{session_num}(~isnan(sensory_enrich_idx{session_num}));
    x = x(~isnan(AUC_data)); % Remove nans from AUC data
    
    [~,temp] = corrcoef(x,y);
    p_value(session_num) = temp(1,2);
    
    % Using simultaneous functional bound
    [fit_curve,gof] = fit(x',y','poly1');
    confi_interval = predint(fit_curve,x,0.95,'functional','on');
    [~,x_idx] = sort(x);
    p=patch([x(x_idx) fliplr(x(x_idx))], [confi_interval(x_idx,1)'  fliplr(confi_interval(x_idx,2)')],[0.8 0.8 0.8],'LineStyle','none');
    p.FaceColor = sessions_color_map(session_num,:);
    p.EdgeColor = 'none';
    p.FaceAlpha = 0.4;
    hold on;
    p1 = plot(x,(x*fit_curve.p1+fit_curve.p2));
    hold on;
    p1(1).Color = sessions_color_map(session_num,:);
    p1(1).LineWidth = 1;
    p1(1).Marker = 'none';
    box(gca,'off')
    legend off
    xlabel('Enrichment index')
    ylabel('Decoding accuracy')
    xlim(x_lim)
    ylim(y_lim)
    yticks([0 0.4 0.8])
    clear x y
end

%% Plot delay
subplot(1,2,2)
for session_num =[1 3]
    
    AUC_data = data.AUC.choice_select_concat_delay1_left_choice_auc{session_num};
    y = AUC_data(~isnan(AUC_data));
    
    % Ensure that the length of AUC data and enrichment indices are the same
    assert(length(AUC_data) == length(delay_enrich_idx{session_num}(~isnan(delay_enrich_idx{session_num}))))
    x = delay_enrich_idx{session_num}(~isnan(delay_enrich_idx{session_num}));
    x = x(~isnan(AUC_data));
    
    [~,temp] = corrcoef(x,y);
    p_value(session_num) = temp(1,2);
    
    % Using simultaneous functional bound
    [fit_curve,gof] = fit(x',y','poly1');
    confi_interval = predint(fit_curve,x,0.95,'functional','on');
    [~,x_idx] = sort(x);
    p=patch([x(x_idx) fliplr(x(x_idx))], [confi_interval(x_idx,1)'  fliplr(confi_interval(x_idx,2)')],[0.8 0.8 0.8],'LineStyle','none');
    p.FaceColor = sessions_color_map(session_num,:);
    p.EdgeColor = 'none';
    p.FaceAlpha = 0.4
    hold on;
    p1 = plot(x,(x*fit_curve.p1+fit_curve.p2))
    hold on;
    p1(1).Color = sessions_color_map(session_num,:);
    p1(1).LineWidth = 1;
    p1(1).Marker = 'none';
    box(gca,'off')
    legend off
    title('Choice')
    xlabel('Enrichment index')
    ylabel('Decoding accuracy')
    xlim(x_lim)
    ylim(y_lim)
    yticks([0 0.4 0.8])
    clear x y
end
end