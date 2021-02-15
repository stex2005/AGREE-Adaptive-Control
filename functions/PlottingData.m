function [OutputPlot]=PlottingData(EMGPreProcessed,DataPreProcessed,BetaFunction_duration,Modality)

WindowNumber=length(DataPreProcessed.Biceps)/BetaFunction_duration;

for window_i=1:WindowNumber
    EMGSignal_Biceps(:,window_i)=EMGPreProcessed.Biceps(1+BetaFunction_duration*(window_i-1):BetaFunction_duration*window_i);
    EMGSignal_Triceps(:,window_i)=EMGPreProcessed.Triceps(1+BetaFunction_duration*(window_i-1):BetaFunction_duration*window_i);
    
    MeasuredPosition(:,window_i)=DataPreProcessed.MeasuredPosition(1+BetaFunction_duration*(window_i-1):BetaFunction_duration*window_i);
    DesiredPosition(:,window_i)=DataPreProcessed.DesiredPosition(1+BetaFunction_duration*(window_i-1):BetaFunction_duration*window_i);
    MeasuredVelocity(:,window_i)=DataPreProcessed.MeasuredVelocity(1+BetaFunction_duration*(window_i-1):BetaFunction_duration*window_i);
    Torque(:,window_i)=DataPreProcessed.Torque(1+BetaFunction_duration*(window_i-1):BetaFunction_duration*window_i);

end

% Mean and std value
OutputPlot.Biceps.SignalMean=mean(EMGSignal_Biceps,2);
OutputPlot.Biceps.SignalStd=std(EMGSignal_Biceps')';

OutputPlot.Triceps.SignalMean=mean(EMGSignal_Triceps,2);
OutputPlot.Triceps.SignalStd=std(EMGSignal_Triceps')';

OutputPlot.Encoder.MeasuredPosition_Mean=mean(MeasuredPosition,2);
OutputPlot.Encoder.MeasuredPosition_Std=std(MeasuredPosition')';

OutputPlot.Encoder.DesiredPosition_Mean=mean(DesiredPosition,2);
OutputPlot.Encoder.DesiredPosition_Std=std(DesiredPosition')';

OutputPlot.Encoder.MeasuredVelocity_Mean=mean(MeasuredVelocity,2);
OutputPlot.Encoder.MeasuredVelocity_Std=std(MeasuredVelocity')';

OutputPlot.Encoder.Torque_Mean=mean(Torque,2);
OutputPlot.Encoder.Torque_Std=std(Torque')';

n_outcome=2;
testo=Modality;

figure,
sgtitle(string(testo),'Interpreter','LaTex')

% Biceps
subplot(n_outcome,1,1)
devpb=OutputPlot.Biceps.SignalMean+OutputPlot.Biceps.SignalStd;
devnb=OutputPlot.Biceps.SignalMean-OutputPlot.Biceps.SignalStd;
plot(devpb,'--k')
hold on
plot(devnb,'--k')
x=1:1:BetaFunction_duration;
x=x';
hold on
plot(OutputPlot.Biceps.SignalMean,'b')
title('Biceps - Normalized EMG Signal','Interpreter','LaTex')
ylabel('$EMG$ $Amplitude$','Interpreter','LaTex')
xlabel('$Time$ [ms]','Interpreter','LaTex')
axesH = gca;
axesH.FontSize=10;
axesH.TickLabelInterpreter='LaTex';
ylim([0 1])

% Triceps
subplot(n_outcome,1,2)
devpb=OutputPlot.Triceps.SignalMean+OutputPlot.Triceps.SignalStd;
devnb=OutputPlot.Triceps.SignalMean-OutputPlot.Triceps.SignalStd;
plot(devpb,'--k')
hold on
plot(devnb,'--k')
x=1:1:BetaFunction_duration;
x=x';
hold on
plot(OutputPlot.Triceps.SignalMean,'b')
title('Triceps - Normalized EMG Signal','Interpreter','LaTex')
ylabel('$EMG$ $Amplitude$','Interpreter','LaTex')
xlabel('$Time$ [ms]','Interpreter','LaTex')
ylim([0 1])
axesH = gca;
axesH.FontSize=10;
axesH.TickLabelInterpreter='LaTex';

% % Position
% subplot(n_outcome,1,3)
% devpb=OutputPlot.Encoder.MeasuredPosition_Mean+OutputPlot.Encoder.MeasuredPosition_Std;
% devnb=OutputPlot.Encoder.MeasuredPosition_Mean-OutputPlot.Encoder.MeasuredPosition_Std;
% plot(devpb,'--k')
% hold on
% plot(devnb,'--k')
% x=1:1:BetaFunction_duration;
% x=x';
% hold on
% plot(OutputPlot.Encoder.MeasuredPosition_Mean,'b')
% hold on
% plot(OutputPlot.Encoder.DesiredPosition_Mean,'r')
% title('Kinematics - Angular position','Interpreter','LaTex')
% ylabel('$Angular$ $position$ [rad]','Interpreter','LaTex')
% xlabel('$Time$ [ms]','Interpreter','LaTex')
% axesH = gca;
% axesH.FontSize=10;
% axesH.TickLabelInterpreter='LaTex';
% 
% % Torque
% subplot(n_outcome,1,4)
% devpb=OutputPlot.Encoder.Torque_Mean+OutputPlot.Encoder.Torque_Std;
% devnb=OutputPlot.Encoder.Torque_Mean-OutputPlot.Encoder.Torque_Std;
% plot(devpb,'--k')
% hold on
% plot(devnb,'--k')
% x=1:1:BetaFunction_duration;
% x=x';
% hold on
% plot(OutputPlot.Encoder.Torque_Mean,'b')
% title('Torque','Interpreter','LaTex')
% ylabel('$Torque$ [mNm]','Interpreter','LaTex')
% xlabel('$Time$ [ms]','Interpreter','LaTex')
% ylim([-4000 2500])
% axesH = gca;
% axesH.FontSize=10;
% axesH.TickLabelInterpreter='LaTex';



end