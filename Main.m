clear all
close all
clc

%% Dependencies
%Adding functions written for this simulation
curr_path = pwd;
functions_lib_path = [curr_path '/functions'];
addpath(functions_lib_path);
functions_lib_path = [curr_path '/onedof_log'];
addpath(functions_lib_path);

%% Data loading

% User: Vale
% load User1-VL.mat

% User: Marta
% load User3-MG.mat

% User: India
load User4-IC.mat

figure
sgtitle('Raw data','Interpreter','LaTex')


subplot(4,2,1),plot(User.Passive.Biceps)
hold on,plot(User.Passive.Triceps)
hold on,plot(-User.Passive.MeasuredPosition)
ylim([-2 4])
title('Passive','Interpreter','LaTex')
legend('Biceps','Triceps','Measured position','Interpreter','LaTex')

subplot(4,2,2),plot(User.AAN.Biceps)
hold on,plot(User.AAN.Triceps)
hold on,plot(-User.AAN.MeasuredPosition)
ylim([-2 4])
title('Partially-assistive (k=10, d=2)','Interpreter','LaTex')
legend('Biceps','Triceps','Measured position','Interpreter','LaTex')

subplot(4,2,3),plot(User.AAN2.Biceps)
hold on,plot(User.AAN2.Triceps)
hold on,plot(-User.AAN2.MeasuredPosition)
ylim([-2 4])
title('Partially-assistive (k=5, d=1)','Interpreter','LaTex')
legend('Biceps','Triceps','Measured position','Interpreter','LaTex')

subplot(4,2,4),plot(User.Antig.Biceps)
hold on,plot(User.Antig.Triceps)
hold on,plot(-User.Antig.MeasuredPosition)
ylim([-2 4])
title('Anti-g','Interpreter','LaTex')
legend('Biceps','Triceps','Measured position','Interpreter','LaTex')

subplot(4,2,5),plot(User.Transparent.Biceps)
hold on,plot(User.Transparent.Triceps)
hold on,plot(-User.Transparent.MeasuredPosition)
ylim([-2 4])
title('Transparent','Interpreter','LaTex')
legend('Biceps','Triceps','Measured position','Interpreter','LaTex')

subplot(4,2,6),plot(User.Resistive.Biceps)
hold on,plot(User.Resistive.Triceps)
hold on,plot(-User.Resistive.MeasuredPosition)
ylim([-2 4])
title('Resistive','Interpreter','LaTex')
legend('Biceps','Triceps','Measured position','Interpreter','LaTex')

subplot(4,2,7),plot(User.Challenging.Biceps)
hold on,plot(User.Challenging.Triceps)
hold on,plot(-User.Challenging.MeasuredPosition)
ylim([-2 4])
title('Challenging','Interpreter','LaTex')
legend('Biceps','Triceps','Measured position','Interpreter','LaTex')


%% Data pre processing
n_repetition=15;
first_repetition=4;
last_repetition=14;
BetaFunction_duration=8000;

Frequency_low=4;
Frequency_high=10;
Fs=1000;

[EMGPreProcessed,DataPreProcessed]=DataPreProcessing(User,first_repetition,last_repetition,BetaFunction_duration,Frequency_low,Frequency_high,Fs);

figure
sgtitle('Pre Processed data','Interpreter','LaTex')
subplot(3,2,1),plot(EMGPreProcessed.Passive.Biceps)
hold on,plot(EMGPreProcessed.Passive.Triceps)
hold on,plot(-User.Passive.MeasuredPosition)
ylim([-2 1])
title('Passive','Interpreter','LaTex')
legend('Biceps','Triceps','Measured position','Interpreter','LaTex')

subplot(3,2,2),plot(EMGPreProcessed.AAN.Biceps)
hold on,plot(EMGPreProcessed.AAN.Triceps)
hold on,plot(-User.AAN.MeasuredPosition)
ylim([-2 1])
title('Partially-assistive','Interpreter','LaTex')
legend('Biceps','Triceps','Measured position','Interpreter','LaTex')

subplot(3,2,3),plot(EMGPreProcessed.Antig.Biceps)
hold on,plot(EMGPreProcessed.Antig.Triceps)
hold on,plot(-User.Antig.MeasuredPosition)
ylim([-2 1])
title('Anti-g','Interpreter','LaTex')
legend('Biceps','Triceps','Measured position','Interpreter','LaTex')

subplot(3,2,4),plot(EMGPreProcessed.Transparent.Biceps)
hold on,plot(EMGPreProcessed.Transparent.Triceps)
hold on,plot(-User.Transparent.MeasuredPosition)
ylim([-2 1])
title('Transparent','Interpreter','LaTex')
legend('Biceps','Triceps','Measured position','Interpreter','LaTex')

subplot(3,2,5),plot(EMGPreProcessed.Resistive.Biceps)
hold on,plot(EMGPreProcessed.Resistive.Triceps)
hold on,plot(-User.Resistive.MeasuredPosition)
ylim([-2 1])
title('Resistive','Interpreter','LaTex')
legend('Biceps','Triceps','Measured position','Interpreter','LaTex')

subplot(3,2,6),plot(EMGPreProcessed.Challenging.Biceps)
hold on,plot(EMGPreProcessed.Challenging.Triceps)
hold on,plot(-User.Challenging.MeasuredPosition)
ylim([-2 1])
title('Challenging','Interpreter','LaTex')
legend('Biceps','Triceps','Measured position','Interpreter','LaTex')

%% Outcome measures computation
[OutcomeMeasure_AllRepetitions.Passive,OutcomeMeasure_Tot{1}]=OutcomeMeasureComputation(EMGPreProcessed.Passive,DataPreProcessed.Passive,BetaFunction_duration);
[OutcomeMeasure_AllRepetitions.AAN,OutcomeMeasure_Tot{2}]=OutcomeMeasureComputation(EMGPreProcessed.AAN,DataPreProcessed.AAN,BetaFunction_duration);
[OutcomeMeasure_AllRepetitions.Antig,OutcomeMeasure_Tot{3}]=OutcomeMeasureComputation(EMGPreProcessed.Antig,DataPreProcessed.Antig,BetaFunction_duration);
[OutcomeMeasure_AllRepetitions.Transparent,OutcomeMeasure_Tot{4}]=OutcomeMeasureComputation(EMGPreProcessed.Transparent,DataPreProcessed.Transparent,BetaFunction_duration);
[OutcomeMeasure_AllRepetitions.Resistive,OutcomeMeasure_Tot{5}]=OutcomeMeasureComputation(EMGPreProcessed.Resistive,DataPreProcessed.Resistive,BetaFunction_duration);
[OutcomeMeasure_AllRepetitions.Challenging,OutcomeMeasure_Tot{6}]=OutcomeMeasureComputation(EMGPreProcessed.Challenging,DataPreProcessed.Challenging,BetaFunction_duration);

%% Plot
Modality='Passive';
OutputPlot.Passive=PlottingData(EMGPreProcessed.Passive,DataPreProcessed.Passive,BetaFunction_duration,Modality);
Modality='Assistive';
OutputPlot.AAN=PlottingData(EMGPreProcessed.AAN,DataPreProcessed.AAN,BetaFunction_duration,Modality);
Modality='Anti-gravity';
OutputPlot.Antig=PlottingData(EMGPreProcessed.Antig,DataPreProcessed.Antig,BetaFunction_duration,Modality);
Modality='Transparent';
OutputPlot.Transparent=PlottingData(EMGPreProcessed.Transparent,DataPreProcessed.Transparent,BetaFunction_duration,Modality);
Modality='Resistive';
OutputPlot.Resistive=PlottingData(EMGPreProcessed.Resistive,DataPreProcessed.Resistive,BetaFunction_duration,Modality);
Modality='Challenging';
OutputPlot.Challenging=PlottingData(EMGPreProcessed.Challenging,DataPreProcessed.Challenging,BetaFunction_duration,Modality);




