clear all 
close all
clc


%% Call Data %%
mat = dir('*.mat'); %1-14 are old participants, 16-31 are young participants 


for i = 1:length(mat);
    cont = load(mat(i).name);
    avg_AnkM_AP{i} = cont.avg_AnkM_AP;
    avg_KneeM_AP{i} = cont.avg_KneeM_AP;
    avg_HipM_AP{i} = cont.avg_HipM_AP;
end


Ank_AP = cell2mat(avg_AnkM_AP);
Knee_AP = cell2mat(avg_KneeM_AP);
Hip_AP = cell2mat(avg_HipM_AP);


y_AnkAP = mean(Ank_AP(:,15:end),2);
o_AnkAP = mean(Ank_AP(:,1:14),2);



y_KneeAP = mean(Knee_AP(:,15:end),2);
o_KneeAP = mean(Knee_AP(:,1:14),2);


y_HipAP = mean(Hip_AP(:,15:end),2);
o_HipAP = mean(Hip_AP(:,1:14),2);

figure(1)
plot(y_AnkAP, "blue")
hold on
plot(o_AnkAP, "red")
hold off

figure(2)
plot(y_KneeAP, "blue")
hold on
plot(o_KneeAP, "red")
hold off


figure(3)
plot(y_HipAP, "blue")
hold on
plot(o_HipAP, "red")
hold off



figure(4)
subplot(2,1,1)
plot(Ank_AP(:,15:end))
subplot(2,1,2)
plot(Ank_AP(:,1:14))



%% Joint ML Moments %% 


for i = 1:length(mat);
    cont = load(mat(i).name);
    avg_AnkM_ML{i} = cont.avg_AnkM_ML;
    avg_KneeM_ML{i} = cont.avg_KneeM_ML;
    avg_HipM_ML{i} = cont.avg_HipM_ML;
end


Ank_ML = cell2mat(avg_AnkM_ML);
Knee_ML = cell2mat(avg_KneeM_ML);
Hip_ML = cell2mat(avg_HipM_ML);


y_AnkML = mean(Ank_ML(:,15:end),2);
o_AnkML = mean(Ank_ML(:,1:14),2);



y_KneeML = mean(Knee_ML(:,15:end),2);
o_KneeML = mean(Knee_ML(:,1:14),2);


y_HipML = mean(Hip_ML(:,15:end),2);
o_HipML = mean(Hip_ML(:,1:14),2);

figure(1)
plot(y_AnkML, "blue")
hold on
plot(o_AnkML, "red")
hold off

figure(2)
plot(y_KneeML, "blue")
hold on
plot(o_KneeML, "red")
hold off


figure(3)
plot(y_HipML, "blue")
hold on
plot(o_HipML, "red")
hold off





figure(4)
subplot(2,1,1)
plot(Hip_ML(:,15:end))
subplot(2,1,2)
plot(Hip_ML(:,1:14))


%% Joint Powers Plotting %% 


for i = 1:length(mat);
    cont = load(mat(i).name);
    avg_AnkP_AP{i} = cont.avg_AnkP_AP;
    avg_KneeP_AP{i} = cont.avg_KneeP_AP;
    avg_HipP_AP{i} = cont.avg_HipP_AP;
end


Ank_Pwr_AP = cell2mat(avg_AnkP_AP);
Knee_Pwr_AP = cell2mat(avg_KneeP_AP);
Hip_Pwr_AP = cell2mat(avg_HipP_AP);


y_AnkPwr_AP = mean(Ank_Pwr_AP(:,15:end),2);
o_AnkPwr_AP = mean(Ank_Pwr_AP(:,1:14),2);



y_KneePwr_AP = mean(Knee_Pwr_AP(:,15:end),2);
o_KneePwr_AP = mean(Knee_Pwr_AP(:,1:14),2);


y_HipPwr_AP = mean(Hip_Pwr_AP(:,15:end),2);
o_HipPwr_AP = mean(Hip_Pwr_AP(:,1:14),2);


figure(1)
plot(y_AnkPwr_AP)
hold on
plot(o_AnkPwr_AP)
hold off

figure(2)
plot(y_KneePwr_AP)
hold on
plot(o_KneePwr_AP)
hold off


figure(3)
plot(y_HipPwr_AP)
hold on
plot(o_HipPwr_AP)
hold off


%% MOS plotting %%

for i = 1:length(mat);
    cont = load(mat(i).name);
    MOS_AntP{i} = cont.MOS_AntP;
    MOS_MedL{i} = cont.MOS_MedL;
end




MOS_AP = cell2mat(MOS_AntP);
MOS_ML = cell2mat(MOS_MedL);

MOS_AP_Y = MOS_AP(:,15:end)';
MOS_ML_Y = MOS_ML(:,15:end)';

MOS_AP_O = MOS_AP(:,1:14)';
MOS_ML_O = MOS_ML(:,1:14)';


MOSAP_old = mean(MOS_AP(:,1:14),2);
MOSAP_yng = mean(MOS_AP(:,15:end),2);

MOSML_old = mean(MOS_ML(:,1:14),2);
MOSML_yng = mean(MOS_ML(1,15:end),2);


avg_MOSAP_Y = mean(MOS_AP_Y,2);
avg_MOSAP_O = mean(MOS_AP_O,2);

avg_MOSML_Y = mean(MOS_ML_Y,2);
avg_MOSML_O = mean(MOS_ML_O,2);


plot1 = MOS_ML_O';
plot2 = MOS_ML_Y';

figure(1)
subplot(2,1,1)
plot(plot1)
subplot(2,1,2)
plot(plot2)

figure(2)
plot(mean(plot1,2), "red", 'LineWidth',2)
hold on
plot(mean(plot2,2),"blue")
hold off 


%% SPM plotting for MOS %%

% Generate some example data (replace this with your actual data)
data = plot1;

% Calculate mean and standard deviation
mean_data = mean(data, 2);
std_data = std(data, 0, 2);



% Calculate mean and standard deviation for the second variable
mean_data2 = mean(plot2, 2);
std_data2 = std(plot2, 0, 2);



% Define x values (assuming you're plotting against the row indices)
x_values = 1:size(data, 1);

% Plot the mean line for the first variable
plot(x_values, mean_data, 'b-', 'LineWidth', 2);
hold on;

% Plot the shaded error region for the first variable
fill([x_values, fliplr(x_values)], [mean_data' + std_data', fliplr(mean_data' - std_data')], 'b', 'FaceAlpha', 0.3, 'EdgeAlpha', 0);

% Plot the mean line for the second variable
plot(x_values, mean_data2, 'r-', 'LineWidth', 2);

% Plot the shaded error region for the second variable
fill([x_values, fliplr(x_values)], [mean_data2' + std_data2', fliplr(mean_data2' - std_data2')], 'r', 'FaceAlpha', 0.3, 'EdgeAlpha', 0);

% Add labels and legend
xlabel('% Single Step');
ylabel('Mediolateral Margin of Stability (m)');
title('Line Graph with Shaded Standard Deviation Error');

% Create legend entries
legend_entries = {'Old Adults Mean', 'Old Adults Standard Deviation', ...
                  'Young Adults Mean', 'Young Adults Standard Deviation', ...
                  'Statistical Significance'};

% Plot legend with legend entries
legend(legend_entries);

% Add rectangle at the bottom for statistical significance
sig_start = 53;
sig_end = 99;
sig_level = 0.05; % significance level (adjust as needed)

if sig_start >= 1 && sig_end <= length(x_values)
    sig_interval = sig_start:sig_end;
    sig_mask = mean_data(sig_interval) + std_data(sig_interval) > mean_data2(sig_interval) - std_data2(sig_interval);
    sig_x = x_values(sig_interval);
    sig_y = 0; % Position of the rectangle at the bottom
    sig_x_start = sig_x(1);
    sig_x_end = sig_x(end);
    rectangle('Position', [sig_x_start, sig_y, sig_x_end - sig_x_start, 0.005], 'FaceColor', 'k', 'EdgeColor', 'none');
    
    % Dummy plot for the legend (outside the visible range)
    plot(NaN, NaN, 's', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerSize', 10); % Square marker for legend
end


% Adjust the y-axis scale
ylim([0, 0.2]); % Set the y-axis limits from -3 to 3
% Adjust the x-axis scale
xlim([0, 101]); % Set the y-axis limits from -3 to 3

% Optionally, you can add grid
%grid on;

% Optionally, customize plot appearance further as needed

% Turn hold off to prevent further plotting on this figure
hold off;
%% SPM MOS ML%% 




spm = spm1d.stats.ttest2(MOS_ML_O, MOS_ML_Y)
spmi = spm.inference(0.05, 'two_tailed', true, 'interp',true)
disp(spmi)


%(2) Plot:
close all
spmi.plot()
spmi.plot_threshold_label();
spmi.plot_p_values();


%% SPM MOS AP%% 
spm = spm1d.stats.ttest2(MOS_AP_O, MOS_AP_Y)
spmi = spm.inference(0.05, 'two_tailed', true, 'interp',true)
disp(spmi)


%(2) Plot:
close all
spmi.plot()
spmi.plot_threshold_label();
spmi.plot_p_values();



%% SPM Moments ML Ank %% 


figure(1)
plot(y_AnkML, "blue")
hold on
plot(o_AnkML, "red")
hold off

O_AnkML = Ank_ML(:,1:14)';
Y_AnkML = Ank_ML(:,15:end)';

spm = spm1d.stats.ttest2(O_AnkML, Y_AnkML)
spmi = spm.inference(0.05, 'two_tailed', true, 'interp',true)
disp(spmi)


%(2) Plot:
close all
spmi.plot()
spmi.plot_threshold_label();
spmi.plot_p_values();

%% SPM Hip Mom ML %% 
figure(1)
plot(y_HipML, "blue")
hold on
plot(o_HipML, "red")
hold off


O_HipML = Hip_ML(:,1:14)';
Y_HipML = Hip_ML(:,15:end)';

spm = spm1d.stats.ttest2(O_HipML, Y_HipML)
spmi = spm.inference(0.05, 'two_tailed', true, 'interp',true)
disp(spmi)


%(2) Plot:
close all
spmi.plot()
spmi.plot_threshold_label();
spmi.plot_p_values();



%% SPM Hip Mom AP %% 
figure(1)
plot(y_HipAP, "blue")
hold on
plot(o_HipAP, "red")
hold off


O_HipAP = Hip_AP(:,1:14)';
Y_HipAP = Hip_AP(:,15:end)';

spm = spm1d.stats.ttest2(O_HipAP, Y_HipAP)
spmi = spm.inference(0.05, 'two_tailed', true, 'interp',true)
disp(spmi)


%(2) Plot:
close all
spmi.plot()
spmi.plot_threshold_label();
spmi.plot_p_values();



%% SPM Ankle Moment AP %% 
figure(1)
plot(y_AnkAP, "blue")
hold on
plot(o_AnkAP, "red")
hold off


O_AnkAP = Ank_AP(:,1:14)';
Y_AnkAP = Ank_AP(:,15:end)';

spm = spm1d.stats.ttest2(O_AnkAP, Y_AnkAP)
spmi = spm.inference(0.05, 'two_tailed', true, 'interp',true)
disp(spmi)


%(2) Plot:
close all
spmi.plot()
spmi.plot_threshold_label();
spmi.plot_p_values();


%% SPM Knee Moment AP %% 

figure(1)
plot(y_KneeAP, "blue")
hold on
plot(o_KneeAP, "red")
hold off


O_KneeAP = Knee_AP(:,1:14)';
Y_KneeAP = Knee_AP(:,15:end)';

spm = spm1d.stats.ttest2(O_KneeAP, Y_KneeAP)
spmi = spm.inference(0.05, 'two_tailed', true, 'interp',true)
disp(spmi)


%(2) Plot:
close all
spmi.plot()
spmi.plot_threshold_label();
spmi.plot_p_values();


%% Cross Correlation Plot Moments against MOS in  AP and ML%% 

hip_AP_yng = Hip_AP(1:60,15:end);
hip_AP_yng = interpft(hip_AP_yng,101);


hip_AP_old = Hip_AP(1:60,1:14);
hip_AP_old = interpft(hip_AP_old,101);

mos_ap_old = (MOS_AP_O')*-1;
mos_ap_yng = (MOS_AP_Y')*-1;



hip_ML_yng = Hip_ML(1:60,15:end);
hip_ML_yng = interpft(hip_ML_yng,101);


hip_ML_old = Hip_ML(1:60,1:14);
hip_ML_old = interpft(hip_ML_old,101);

mos_ml_old = (MOS_ML_O');
mos_ml_yng = (MOS_ML_Y');



ank_ML_yng = Ank_ML(1:60,15:end);
ank_ML_yng = interpft(ank_ML_yng,101);


ank_ML_old = Ank_ML(1:60,1:14);
ank_ML_old = interpft(ank_ML_old,101);



ank_AP_yng = Ank_AP(1:60,15:end);
ank_AP_yng = interpft(ank_AP_yng,101);


ank_AP_old = Ank_AP(1:60,1:14);
ank_AP_old = interpft(ank_AP_old,101);




knee_AP_yng = Knee_AP(1:60,15:end);
knee_AP_yng = interpft(knee_AP_yng,101);


knee_AP_old = Knee_AP(1:60,1:14);
knee_AP_old = interpft(knee_AP_old,101);


Syn_yng = [hip_AP_yng + knee_AP_yng + ank_AP_yng];
Syn_old = [hip_AP_old + knee_AP_old + ank_AP_old];


for i = 1:size(hip_AP_yng,2);
        [r_yng(:,i), lags_yng(:,i)] = crosscorr(hip_AP_yng(:,i), mos_ap_yng(:,i), NumLags=50);
        [r_MLyng(:,i), lags_MLyng(:,i)] = crosscorr(hip_ML_yng(:,i), mos_ml_yng(:,i), NumLags=50);
        [r_Ank_MLyng(:,i), lags_Ank_MLyng(:,i)] = crosscorr(ank_ML_yng(:,i), mos_ml_yng(:,i), NumLags=50);
        [r_Ank_APyng(:,i), lags_Ank_APyng(:,i)] = crosscorr(ank_AP_yng(:,i), mos_ap_yng(:,i), NumLags=50);
        [r_Synyng(:,i), lags_Synyng(:,i)] = crosscorr(Syn_yng(:,i), mos_ap_yng(:,i), NumLags=50);
        [r_Knee_APyng(:,i), lags_Knee_APyng(:,i)] = crosscorr(knee_AP_yng(:,i), mos_ap_yng(:,i), NumLags=50);
end 


for i = 1:size(hip_AP_old,2)
     [r_old(:,i), lags_old(:,i)] = crosscorr(hip_AP_old(:,i), mos_ap_old(:,i), NumLags=50);
     [r_MLold(:,i), lags_MLold(:,i)] = crosscorr(hip_ML_old(:,i), mos_ml_old(:,i), NumLags=50);
     [r_Ank_MLold(:,i), lags_Ank_MLold(:,i)] = crosscorr(ank_ML_old(:,i), mos_ml_old(:,i), NumLags=50);
     [r_Ank_APold(:,i), lags_Ank_APold(:,i)] = crosscorr(ank_AP_old(:,i), mos_ap_old(:,i), NumLags=50);
     [r_Synold(:,i), lags_Synold(:,i)] = crosscorr(Syn_old(:,i), mos_ap_old(:,i), NumLags = 50);
     [r_Knee_APold(:,i), lags_Knee_APold(:,i)] = crosscorr(knee_AP_old(:,i), mos_ap_old(:,i), NumLags=50);

end 




for i = 1:size(r_yng,2);
    % Create a new figure
    figure;
    
    % Plot something (optional)
    stem(lags_MLyng(:,i),r_MLyng(:,i));
end

% Find zero values in the matrix
%zero_HipAP_indices = lags_yng == 0;
%zero_HipML_indices = lags_MLyng == 0;

% Display the indices of zero values
%[row_HipAP_indices, column_HipAP_indices] = find(zero_HipAP_indices);
%[row_HipML_indices, column_HipML_indices] = find(zero_HipML_indices);

%row_HipAP_indices= row_HipAP_indices';
%row_HipML_indices= row_HipML_indices';


R_HipAP_yng = r_yng(51,:);
avg_RHipAP_y = mean(R_HipAP_yng)
R_HipAP_old = r_old(51,:);
avg_RHipAP_o = mean(R_HipAP_old)

R_HipML_yng = r_MLyng(51,:);
avg_RHipML_y = mean(R_HipML_yng)
R_HipML_old = r_MLold(51,:);
avg_RHipML_o = mean(R_HipML_old)


R_AnkAP_yng = r_Ank_APyng(51,:);
avg_RAnkAP_y = mean(R_AnkAP_yng)
R_AnkAP_old = r_Ank_APold(51,:);
avg_RAnkAP_o = mean(R_AnkAP_old)


R_AnkML_yng = r_Ank_MLyng(51,:);
avg_RAnkML_y = mean(R_AnkML_yng)
R_AnkML_old = r_Ank_MLold(51,:);
avg_RAnkML_o = mean(R_AnkML_old)


R_Syn_yng = r_Synyng(51,:);
avg_RSyn_yng = mean(R_Syn_yng)
R_Syn_old = r_Synold(51,:);
avg_RSyn_old = mean(R_Syn_old)


R_KneeAP_yng = r_Knee_APyng(51,:);
avg_RKneeAP_y = mean(R_KneeAP_yng)
R_KneeAP_old = r_Knee_APold(51,:);
avg_RKneeAP_o = mean(R_KneeAP_old)


