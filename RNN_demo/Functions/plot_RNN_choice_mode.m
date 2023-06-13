%% Plot choice mode with distractor trials
function plot_RNN_choice_mode(RNN_data,example_RNN)
ppc_left_choice_mode_std_norm = RNN_data.PPC.choice_mode.left_choice_mode_std_norm{example_RNN};
ppc_right_choice_mode_std_norm = RNN_data.PPC.choice_mode.right_choice_mode_std_norm{example_RNN};
ppc_L_pos_choice_mode_std_norm = RNN_data.PPC.choice_mode.L_pos_choice_mode_std_norm{example_RNN};

vs1_left_choice_mode_std_norm  = RNN_data.vS1.choice_mode.left_choice_mode_std_norm{example_RNN};
vs1_right_choice_mode_std_norm = RNN_data.vS1.choice_mode.right_choice_mode_std_norm{example_RNN};
vs1_L_pos_choice_mode_std_norm = RNN_data.vS1.choice_mode.L_pos_choice_mode_std_norm{example_RNN};

PPC_L_pos_mid_point = RNN_data.PPC.choice_mode.L_pos_mid_point{example_RNN};
vS1_L_pos_mid_point = RNN_data.vS1.choice_mode.L_pos_mid_point{example_RNN};

fs_image = 93.5;

%% PPC
x_lim = [1 round(fs_image*3.5)];
left_choice_mode_std_norm = ppc_left_choice_mode_std_norm;
right_choice_mode_std_norm = ppc_right_choice_mode_std_norm;
L_pos_choice_mode_std_norm = ppc_L_pos_choice_mode_std_norm;
rand_idx = randsample(1:length(PPC_L_pos_mid_point),40); % Randomly select 40 trials out of 100

figure('Position',[200,100,300,150],'Color','white','DefaultAxesFontSize',14);
subplot(1,2,1); hold on;
plot(nanmean(left_choice_mode_std_norm,1),'Color',[0.00,0.45,0.74],'LineWidth',1);
plot(nanmean(right_choice_mode_std_norm,1),'Color',[0.64,0.08,0.18],'LineWidth',1);
for j = rand_idx
    if PPC_L_pos_mid_point(j) < 0.5
        plot(L_pos_choice_mode_std_norm((j),:),'Color',[0.00,0.00,0.00],'LineWidth',1);
    else
        plot(L_pos_choice_mode_std_norm((j),:),':','Color',[0.00,0.00,0.00],'LineWidth',1);
    end
end
line([0.5.*fs_image,0.5.*fs_image],[-0.3 1],'Color',[0.25,0.25,0.25]);
line([1.5.*fs_image,1.5.*fs_image],[-0.3 1],'Color',[0.25,0.25,0.25]);
ylabel('Choice activity')
xlabel('Time to go cue')
xticks([])
yticks([0 1])
ylim([-0.3 1])
xlim(x_lim)
title('+PPC')

%% vS1
left_choice_mode_std_norm = vs1_left_choice_mode_std_norm;
right_choice_mode_std_norm = vs1_right_choice_mode_std_norm;
L_pos_choice_mode_std_norm = vs1_L_pos_choice_mode_std_norm;
rand_idx = randsample(1:length(vS1_L_pos_mid_point),40); % Randomly select 40 trials out of 100
subplot(1,2,2); hold on;
plot(nanmean(left_choice_mode_std_norm,1),'Color',[0.00,0.45,0.74],'LineWidth',1);
plot(nanmean(right_choice_mode_std_norm,1),'Color',[0.64,0.08,0.18],'LineWidth',1);
for j = rand_idx
    
    if vS1_L_pos_mid_point(j) < 0.5
        plot(L_pos_choice_mode_std_norm((j),:),'Color',[0.00,0.00,0.00],'LineWidth',1);
    else
        plot(L_pos_choice_mode_std_norm((j),:),':','Color',[0.00,0.00,0.00],'LineWidth',1);
    end
end
line([0.5.*fs_image,0.5.*fs_image],[-0.3 1],'Color',[0.25,0.25,0.25]);
line([1.5.*fs_image,1.5.*fs_image],[-0.3 1],'Color',[0.25,0.25,0.25]);
xlabel('Time to go cue')
xticks([])
yticks([0 1])
xlim(x_lim)
ylim([-0.3 1])
title('+vS1')
end