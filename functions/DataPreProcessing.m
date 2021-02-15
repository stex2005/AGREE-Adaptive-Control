function [EMGPreProcessed,UserPreProcessed]=DataPreProcessing(User,first_repetition,last_repetition,BetaFunction_duration,Frequency_low,Frequency_high,Fs)

%% Select desired repetitions
sync_Modality_start=1+BetaFunction_duration*first_repetition;
sync_Modality_end=BetaFunction_duration*last_repetition;

%% Encoder data cut
UserPreProcessed.Passive.MeasuredPosition=User.Passive.MeasuredPosition(sync_Modality_start:sync_Modality_end,:);
UserPreProcessed.Passive.DesiredPosition=User.Passive.DesiredPosition(sync_Modality_start:sync_Modality_end,:);
UserPreProcessed.Passive.MeasuredVelocity=User.Passive.MeasuredVelocity(sync_Modality_start:sync_Modality_end,:);
UserPreProcessed.Passive.Torque=User.Passive.Torque(sync_Modality_start:sync_Modality_end,:);

UserPreProcessed.AAN.MeasuredPosition=User.AAN.MeasuredPosition(sync_Modality_start:sync_Modality_end,:);
UserPreProcessed.AAN.DesiredPosition=User.AAN.DesiredPosition(sync_Modality_start:sync_Modality_end,:);
UserPreProcessed.AAN.MeasuredVelocity=User.AAN.MeasuredVelocity(sync_Modality_start:sync_Modality_end,:);
UserPreProcessed.AAN.Torque=User.AAN.Torque(sync_Modality_start:sync_Modality_end,:);

UserPreProcessed.Antig.MeasuredPosition=User.Antig.MeasuredPosition(sync_Modality_start:sync_Modality_end,:);
UserPreProcessed.Antig.DesiredPosition=User.Antig.DesiredPosition(sync_Modality_start:sync_Modality_end,:);
UserPreProcessed.Antig.MeasuredVelocity=User.Antig.MeasuredVelocity(sync_Modality_start:sync_Modality_end,:);
UserPreProcessed.Antig.Torque=User.Antig.Torque(sync_Modality_start:sync_Modality_end,:);

UserPreProcessed.Transparent.MeasuredPosition=User.Transparent.MeasuredPosition(sync_Modality_start:sync_Modality_end,:);
UserPreProcessed.Transparent.DesiredPosition=User.Transparent.DesiredPosition(sync_Modality_start:sync_Modality_end,:);
UserPreProcessed.Transparent.MeasuredVelocity=User.Transparent.MeasuredVelocity(sync_Modality_start:sync_Modality_end,:);
UserPreProcessed.Transparent.Torque=User.Transparent.Torque(sync_Modality_start:sync_Modality_end,:);

UserPreProcessed.Resistive.MeasuredPosition=User.Resistive.MeasuredPosition(sync_Modality_start:sync_Modality_end,:);
UserPreProcessed.Resistive.DesiredPosition=User.Resistive.DesiredPosition(sync_Modality_start:sync_Modality_end,:);
UserPreProcessed.Resistive.MeasuredVelocity=User.Resistive.MeasuredVelocity(sync_Modality_start:sync_Modality_end,:);
UserPreProcessed.Resistive.Torque=User.Resistive.Torque(sync_Modality_start:sync_Modality_end,:);

UserPreProcessed.Challenging.MeasuredPosition=User.Challenging.MeasuredPosition(sync_Modality_start:sync_Modality_end,:);
UserPreProcessed.Challenging.DesiredPosition=User.Challenging.DesiredPosition(sync_Modality_start:sync_Modality_end,:);
UserPreProcessed.Challenging.MeasuredVelocity=User.Challenging.MeasuredVelocity(sync_Modality_start:sync_Modality_end,:);
UserPreProcessed.Challenging.Torque=User.Challenging.Torque(sync_Modality_start:sync_Modality_end,:);

%% EMG data cut
UserPreProcessed.Passive.Biceps=User.Passive.Biceps(sync_Modality_start:sync_Modality_end,:);
UserPreProcessed.Passive.Triceps=User.Passive.Triceps(sync_Modality_start:sync_Modality_end,:);

UserPreProcessed.AAN.Biceps=User.AAN.Biceps(sync_Modality_start:sync_Modality_end,:);
UserPreProcessed.AAN.Triceps=User.AAN.Triceps(sync_Modality_start:sync_Modality_end,:);

UserPreProcessed.Antig.Biceps=User.Antig.Biceps(sync_Modality_start:sync_Modality_end,:);
UserPreProcessed.Antig.Triceps=User.Antig.Triceps(sync_Modality_start:sync_Modality_end,:);

UserPreProcessed.Transparent.Biceps=User.Transparent.Biceps(sync_Modality_start:sync_Modality_end,:);
UserPreProcessed.Transparent.Triceps=User.Transparent.Triceps(sync_Modality_start:sync_Modality_end,:);

UserPreProcessed.Resistive.Biceps=User.Resistive.Biceps(sync_Modality_start:sync_Modality_end,:);
UserPreProcessed.Resistive.Triceps=User.Resistive.Triceps(sync_Modality_start:sync_Modality_end,:);

UserPreProcessed.Challenging.Biceps=User.Challenging.Biceps(sync_Modality_start:sync_Modality_end,:);
UserPreProcessed.Challenging.Triceps=User.Challenging.Triceps(sync_Modality_start:sync_Modality_end,:);

%% EMG pre-processing
EMGPreProcessed.Passive.Biceps=EMGFilter(User.Passive.Biceps(sync_Modality_start:sync_Modality_end,:),Frequency_low,Frequency_high,Fs);
EMGPreProcessed.Passive.Triceps=EMGFilter(User.Passive.Triceps(sync_Modality_start:sync_Modality_end,:),Frequency_low,Frequency_high,Fs);

EMGPreProcessed.AAN.Biceps=EMGFilter(User.AAN.Biceps(sync_Modality_start:sync_Modality_end,:),Frequency_low,Frequency_high,Fs);
EMGPreProcessed.AAN.Triceps=EMGFilter(User.AAN.Triceps(sync_Modality_start:sync_Modality_end,:),Frequency_low,Frequency_high,Fs);

EMGPreProcessed.Antig.Biceps=EMGFilter(User.Antig.Biceps(sync_Modality_start:sync_Modality_end,:),Frequency_low,Frequency_high,Fs);
EMGPreProcessed.Antig.Triceps=EMGFilter(User.Antig.Triceps(sync_Modality_start:sync_Modality_end,:),Frequency_low,Frequency_high,Fs);

EMGPreProcessed.Transparent.Biceps=EMGFilter(User.Transparent.Biceps(sync_Modality_start:sync_Modality_end,:),Frequency_low,Frequency_high,Fs);
EMGPreProcessed.Transparent.Triceps=EMGFilter(User.Transparent.Triceps(sync_Modality_start:sync_Modality_end,:),Frequency_low,Frequency_high,Fs);

EMGPreProcessed.Resistive.Biceps=EMGFilter(User.Resistive.Biceps(sync_Modality_start:sync_Modality_end,:),Frequency_low,Frequency_high,Fs);
EMGPreProcessed.Resistive.Triceps=EMGFilter(User.Resistive.Triceps(sync_Modality_start:sync_Modality_end,:),Frequency_low,Frequency_high,Fs);

EMGPreProcessed.Challenging.Biceps=EMGFilter(User.Challenging.Biceps(sync_Modality_start:sync_Modality_end,:),Frequency_low,Frequency_high,Fs);
EMGPreProcessed.Challenging.Triceps=EMGFilter(User.Challenging.Triceps(sync_Modality_start:sync_Modality_end,:),Frequency_low,Frequency_high,Fs);

%% Normalization
AllBicepValue=vertcat(EMGPreProcessed.Passive.Biceps,EMGPreProcessed.AAN.Biceps,EMGPreProcessed.Antig.Biceps,EMGPreProcessed.Transparent.Biceps,EMGPreProcessed.Resistive.Biceps,EMGPreProcessed.Challenging.Biceps);
Maximum_biceps=max(AllBicepValue);

AllTricepValue=vertcat(EMGPreProcessed.Passive.Triceps,EMGPreProcessed.AAN.Triceps,EMGPreProcessed.Antig.Triceps,EMGPreProcessed.Transparent.Triceps,EMGPreProcessed.Resistive.Triceps,EMGPreProcessed.Challenging.Triceps);
Maximum_triceps=max(AllTricepValue);

EMGPreProcessed.Passive.Biceps=EMGPreProcessed.Passive.Biceps/Maximum_biceps;
EMGPreProcessed.Passive.Triceps=EMGPreProcessed.Passive.Triceps/Maximum_triceps;

EMGPreProcessed.AAN.Biceps=EMGPreProcessed.AAN.Biceps/Maximum_biceps;
EMGPreProcessed.AAN.Triceps=EMGPreProcessed.AAN.Triceps/Maximum_triceps;

EMGPreProcessed.Antig.Biceps=EMGPreProcessed.Antig.Biceps/Maximum_biceps;
EMGPreProcessed.Antig.Triceps=EMGPreProcessed.Antig.Triceps/Maximum_triceps;

EMGPreProcessed.Transparent.Biceps=EMGPreProcessed.Transparent.Biceps/Maximum_biceps;
EMGPreProcessed.Transparent.Triceps=EMGPreProcessed.Transparent.Triceps/Maximum_triceps;

EMGPreProcessed.Resistive.Biceps=EMGPreProcessed.Resistive.Biceps/Maximum_biceps;
EMGPreProcessed.Resistive.Triceps=EMGPreProcessed.Resistive.Triceps/Maximum_triceps;

EMGPreProcessed.Challenging.Biceps=EMGPreProcessed.Challenging.Biceps/Maximum_biceps;
EMGPreProcessed.Challenging.Triceps=EMGPreProcessed.Challenging.Triceps/Maximum_triceps;

end