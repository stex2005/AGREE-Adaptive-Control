% Computed parameters:
% EMG - TIME DOMAIN
% Mean absolute value
% Integrated EMG
% Variance of EMG
% Root Mean Square
% Waveform Length
% Skewness

function [Output_AllRepetitions, Output_Mean]=OutcomeMeasureComputation(EMGPreProcessed,DataPreProcessed,BetaFunction_duration)

WindowNumber=length(DataPreProcessed.Biceps)/BetaFunction_duration;

%% 1.1) EMG: Output computation for each single repetitions
for window_i=1:WindowNumber
    % Output to be computed on raw EMG
    
    % Mean absolute value
    Output_AllRepetitions.Biceps.mav(window_i)=mean(abs(DataPreProcessed.Biceps(1+BetaFunction_duration*(window_i-1):BetaFunction_duration*window_i)));
    Output_AllRepetitions.Triceps.mav(window_i)=mean(abs(DataPreProcessed.Triceps(1+BetaFunction_duration*(window_i-1):BetaFunction_duration*window_i)));
    
    % Integrated EMG
    Output_AllRepetitions.Biceps.iemg(window_i)=trapz(abs(DataPreProcessed.Biceps(1+BetaFunction_duration*(window_i-1):BetaFunction_duration*window_i)));
    Output_AllRepetitions.Triceps.iemg(window_i)=trapz(abs(DataPreProcessed.Triceps(1+BetaFunction_duration*(window_i-1):BetaFunction_duration*window_i)));
    
    % Variance of EMG
    Output_AllRepetitions.Biceps.var(window_i)=var(DataPreProcessed.Biceps(1+BetaFunction_duration*(window_i-1):BetaFunction_duration*window_i));
    Output_AllRepetitions.Triceps.var(window_i)=var(DataPreProcessed.Triceps(1+BetaFunction_duration*(window_i-1):BetaFunction_duration*window_i));
    
    % Root Mean Square
    Output_AllRepetitions.Biceps.rms(window_i)=rms(DataPreProcessed.Biceps(1+BetaFunction_duration*(window_i-1):BetaFunction_duration*window_i));
    Output_AllRepetitions.Triceps.rms(window_i)=rms(DataPreProcessed.Triceps(1+BetaFunction_duration*(window_i-1):BetaFunction_duration*window_i));
    
    % Waveform Length
    n=length(DataPreProcessed.Biceps(1+BetaFunction_duration*(window_i-1):BetaFunction_duration*window_i));
    wl_b=0;
    wl_t=0;
    for z=2:n
        wl_b=wl_b+abs(DataPreProcessed.Biceps(1+BetaFunction_duration*(window_i-1)+(z-1))-DataPreProcessed.Biceps(1+BetaFunction_duration*(window_i-1)+(z-2)));
        wl_t=wl_b+abs(DataPreProcessed.Triceps(1+BetaFunction_duration*(window_i-1)+(z-1))-DataPreProcessed.Triceps(1+BetaFunction_duration*(window_i-1)+(z-2)));
    end
    Output_AllRepetitions.Biceps.wl(window_i)=wl_b;
    Output_AllRepetitions.Triceps.wl(window_i)=wl_t;
    
end

%% 1.2) EMG: Mean output computation among all repetitions
% So che non ti piace ma ero in zero-sbatta-style

Output_Mean.Biceps_mav(1)=mean(Output_AllRepetitions.Biceps.mav);
Output_Mean.Biceps_mav(2)=std(Output_AllRepetitions.Biceps.mav);
Output_Mean.Triceps_mav(1)=mean(Output_AllRepetitions.Triceps.mav);
Output_Mean.Triceps_mav(2)=std(Output_AllRepetitions.Triceps.mav);

Output_Mean.Biceps_iemg(1)=mean(Output_AllRepetitions.Biceps.iemg);
Output_Mean.Biceps_iemg(2)=std(Output_AllRepetitions.Biceps.iemg);
Output_Mean.Triceps_iemg(1)=mean(Output_AllRepetitions.Triceps.iemg);
Output_Mean.Triceps_iemg(2)=std(Output_AllRepetitions.Triceps.iemg);

Output_Mean.Biceps_var(1)=mean(Output_AllRepetitions.Biceps.var);
Output_Mean.Biceps_var(2)=std(Output_AllRepetitions.Biceps.var);
Output_Mean.Triceps_var(1)=mean(Output_AllRepetitions.Triceps.var);
Output_Mean.Triceps_var(2)=std(Output_AllRepetitions.Triceps.var);

Output_Mean.Biceps_rms(1)=mean(Output_AllRepetitions.Biceps.rms);
Output_Mean.Biceps_rms(2)=std(Output_AllRepetitions.Biceps.rms);
Output_Mean.Triceps_rms(1)=mean(Output_AllRepetitions.Triceps.rms);
Output_Mean.Triceps_rms(2)=std(Output_AllRepetitions.Triceps.rms);

Output_Mean.Biceps_wl(1)=mean(Output_AllRepetitions.Biceps.wl);
Output_Mean.Biceps_wl(2)=std(Output_AllRepetitions.Biceps.wl);
Output_Mean.Triceps_wl(1)=mean(Output_AllRepetitions.Triceps.wl);
Output_Mean.Triceps_wl(2)=std(Output_AllRepetitions.Triceps.wl);

%% 2.1) Kinematics: Output computation for each single repetitions

for window_i=1:WindowNumber
    
    % Integral error along trajectory
    Output_AllRepetitions.Encoder.integral(window_i)=abs(trapz(DataPreProcessed.DesiredPosition(1+BetaFunction_duration*(window_i-1):BetaFunction_duration*window_i))-trapz(DataPreProcessed.MeasuredPosition(1+BetaFunction_duration*(window_i-1):BetaFunction_duration*window_i)));

    % Mean deviation from desired trajectory
    deviation=abs(DataPreProcessed.DesiredPosition(1+BetaFunction_duration*(window_i-1):BetaFunction_duration*window_i)-DataPreProcessed.MeasuredPosition(1+BetaFunction_duration*(window_i-1):BetaFunction_duration*window_i));
    Output_AllRepetitions.Encoder.ydev(window_i)=sum(deviation)/BetaFunction_duration;
    
end
    
%% 2.2) Kinematics: Mean output computation among all repetitions
 
Output_Mean.Encoder_integral(1)=mean(Output_AllRepetitions.Encoder.integral);
Output_Mean.Encoder_integral(2)=std(Output_AllRepetitions.Encoder.integral);

Output_Mean.Encoder_ydev(1)=mean(Output_AllRepetitions.Encoder.ydev);
Output_Mean.Encoder_ydev(2)=std(Output_AllRepetitions.Encoder.ydev);
    

%% 3.1) Dynamics: Output computation for each single repetitions

for window_i=1:WindowNumber
    
    % Torque
    Output_AllRepetitions.Encoder.torque(window_i)=mean(DataPreProcessed.Torque(1+BetaFunction_duration*(window_i-1):BetaFunction_duration*window_i));
    
end
    
%% 3.2) Dynamics: Mean output computation among all repetitions
 
Output_Mean.Encoder_torque(1)=mean(Output_AllRepetitions.Encoder.torque);
Output_Mean.Encoder_torque(2)=std(Output_AllRepetitions.Encoder.torque);

end