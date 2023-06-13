%% plot choice mode
function plot_choice_mode(parameters, retained_elim_data)
% Initialise variables
retained_right_choice_mode = retained_elim_data.choice_mode.retained_right_choice_mode;
lost_right_choice_mode = retained_elim_data.choice_mode.lost_right_choice_mode;
retained_left_choice_mode = retained_elim_data.choice_mode.retained_left_choice_mode;
lost_left_choice_mode = retained_elim_data.choice_mode.lost_left_choice_mode;

cMap = colormap(lines(5)); % Set colormap
fs_image = parameters.fs_image;
x_lim = [round(fs_image*0),round(fs_image*4.0)]; % Stimulus to end of delay
y_lim = [-0.2 1.5];
ITI_epoch = parameters.ITI_epoch;
session_num = 2; % Intermediate session
animal_idx = 1:13;

figure('Position',[200,100,170,150],'Color','white','DefaultAxesFontSize',14);
for hide = 1
    
    concat_a = nan(size(retained_right_choice_mode{3}{session_num},2),13);
    concat_b = nan(size(retained_right_choice_mode{3}{session_num},2),13);
    concat_c = nan(size(retained_right_choice_mode{3}{session_num},2),13);
    concat_d = nan(size(retained_right_choice_mode{3}{session_num},2),13);
    
    all_concat_a{session_num} = nan(size(retained_right_choice_mode{3}{session_num},2),13);
    all_concat_b{session_num} = nan(size(retained_right_choice_mode{3}{session_num},2),13);
    all_concat_c{session_num} = nan(size(retained_right_choice_mode{3}{session_num},2),13);
    all_concat_d{session_num} = nan(size(retained_right_choice_mode{3}{session_num},2),13);
    
    for animal_num = animal_idx%find(abs(cellfun(@isempty,retained_right_choice_mode)-1))
        if ~isempty(retained_right_choice_mode{animal_num})
            if ~isempty(retained_right_choice_mode{animal_num}{session_num})
                for shuffle_num = 1:size(retained_right_choice_mode{animal_num}{session_num},3)
                    % average across trials
                    a = mean(retained_right_choice_mode{animal_num}{session_num}(:,:,shuffle_num),1);
                    b = mean(lost_right_choice_mode{animal_num}{session_num}(:,:,shuffle_num),1);
                    c = mean(retained_left_choice_mode{animal_num}{session_num}(:,:,shuffle_num),1);
                    d = mean(lost_left_choice_mode{animal_num}{session_num}(:,:,shuffle_num),1);
                    
                    % baseline substract ITI_epoch
                    baseline_a = mean(mean(retained_right_choice_mode{animal_num}{session_num}(:,ITI_epoch,shuffle_num),2),1);
                    baseline_b = mean(mean(lost_right_choice_mode{animal_num}{session_num}(:,ITI_epoch,shuffle_num),2),1);
                    baseline_c = mean(mean(retained_left_choice_mode{animal_num}{session_num}(:,ITI_epoch,shuffle_num),2),1);
                    baseline_d = mean(mean(lost_left_choice_mode{animal_num}{session_num}(:,ITI_epoch,shuffle_num),2),1);
                    
                    a = a-baseline_a;
                    b = b-baseline_b;
                    c = c-baseline_c;
                    d = d-baseline_d;
                    
                    % Get selectivity by right-left
                    shuffle_a(:,shuffle_num) = a-c;
                    shuffle_b(:,shuffle_num) = b-d;
                    
                    clear a b c d baseline_a baseline_b  baseline_c baseline_d
                end
                % concat across animals
                concat_a(:,animal_num) = mean(shuffle_a,2);
                concat_b(:,animal_num) = mean(shuffle_b,2);
                
                all_concat_a{session_num}(:,animal_num) = mean(shuffle_a,2);
                all_concat_b{session_num}(:,animal_num) = mean(shuffle_b,2);
            end
        end
    end
    
    hold on;
    % plot retained activity
    plot(nanmedian(concat_a,2),'Color',cMap(2,:),'LineWidth',2);
    % plot eliminated activity
    plot(nanmedian(concat_b,2),'Color',cMap(1,:),'LineWidth',2);
    box(gca,'off')
    xticks([])
    xlim([x_lim])
    ylim([y_lim])
    vline(round(fs_image*1),':k');
    vline(round(fs_image*2),':k');
    ylabel('Choice activity')
    xlabel('Time to go cue')
    title('Choice mode')
end

% Plot unity line figure
session_num = 2;
choice_epoch_idx = round(fs_image*3):round(fs_image*4);
for hide = 1
    figure('Position',[200,100,150,150],'Color','white','DefaultAxesFontSize',14);
    a = nanmean(all_concat_a{session_num}(choice_epoch_idx,:),1);
    b = nanmean(all_concat_b{session_num}(choice_epoch_idx,:),1);
    hold on
    scatter(b,a,'MarkerFaceColor','k','MarkerEdgeColor','none')
    plot_cross_a = b;
    plot_cross_b = a;
    SEM_a = nanstd(plot_cross_a)./sqrt(length(plot_cross_a));
    SEM_b = nanstd(plot_cross_b)./sqrt(length(plot_cross_b));
    line([nanmean(plot_cross_a)-SEM_a nanmean(plot_cross_a)+SEM_a], [nanmean(plot_cross_b) nanmean(plot_cross_b)],'Color',[0.85 0.325 0.098],'LineWidth',2)
    line([nanmean(plot_cross_a) nanmean(plot_cross_a)], [nanmean(plot_cross_b)-SEM_b nanmean(plot_cross_b)+SEM_b],'Color',[0.85 0.325 0.098],'LineWidth',2)
    xlim([-0.5 9])
    ylim([-0.5 9])
    xticks([0 4 8])
    yticks([0 4 8])
    line([-2 9],[-2 9],'Color','k','LineStyle',':')
    ylabel('Activity (retained)')
    xlabel('Activity (eliminated)')
    axis square
end

% Get significance using bootstrap
for shuffle_num = 1:1000
    rand_idx = randsample(length(a),length(a),'true' );
    shuffle_diff(shuffle_num) = nanmean(a(rand_idx) - b(rand_idx));
end
p_value = sum(shuffle_diff<0)/1000;
end