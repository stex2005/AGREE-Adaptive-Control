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

%% Import EMG
[filename, pathname] = uigetfile('*.csv',['Select a file: ']);
if filename == 0
    disp('No file input. Stopping simulation.');
else
    data=import_emg(filename);
end

data=table2array(data);
EMGData=data(2:end,3:4);
EMGCommand=data(2:end,1);
n_muscle=2;

%% Import Encoder
[filename, pathname] = uigetfile('*.csv',['Select a file: ']);
if filename == 0
    disp('No file input. Stopping simulation.');
else
    data=import_onedof_log(filename);
end

data=table2array(data);
EncoderData=data;
%% EMG raw signal plot

figure();
i=0;

i_max = n_muscle;
i=i+1;
subplot(i_max,1,i)
plot(EMGData(:,1))
hold on
plot(EMGCommand)
legend('Biceps','EMG command');

i=i+1;
subplot(i_max,1,i)
plot(EMGData(:,2))
hold on
plot(EMGCommand)
legend('Triceps','EMG command');

%% Encoder plot

figure
plot(EncoderData(:,4))
hold on, plot(EncoderData(:,10))
hold on, plot(EncoderData(:,1))
legend('Measured position','Desired position','Encoder command')


%% Synchronization
% Computation of the difference in time of the onset of the freeze mode
% (Command 6) in the EMG-related data and the Encoder-related data

EMGsync=find(EMGCommand==101); % 6 refers to freeze status
EMGsync=EMGsync(1);

Encodersync=find(EncoderData(:,1)==101); % 6 refers to freeze status
Encodersync=Encodersync(1);

sync_difference=Encodersync-EMGsync;

% Cut of the first n=sync_difference samples in the Encoder-related data
EncoderData(1:sync_difference,:)=[];


MeasuredPosition=EncoderData(:,4);
DesiredPosition=EncoderData(:,10);

figure,
plot(EMGData(:,1))
hold on,plot(EMGData(:,2))
hold on,plot(MeasuredPosition)
hold on,plot(EMGCommand)
hold on,plot(EncoderData(:,1))
legend('Biceps','Triceps','Measured position','EMG Command','Encoder command')


%% Segmentation according to control modality
try
    % 0) MVC
    EMGcut=find(EMGCommand==6);
    User.MVC.Biceps=EMGData(EMGcut(1):EMGcut(end),1);
    User.MVC.Triceps=EMGData(EMGcut(1):EMGcut(end),2);
catch exception
    disp('MVC not found');
end
% 1) Passive
EMGcut=find(EMGCommand==101);
User.Passive.Biceps=EMGData(EMGcut(1):EMGcut(end),1);
User.Passive.Triceps=EMGData(EMGcut(1):EMGcut(end),2);

Encodercut=find(EncoderData(:,1)==101);
User.Passive.MeasuredPosition=EncoderData(Encodercut(1):Encodercut(end),4);
User.Passive.DesiredPosition=EncoderData(Encodercut(1):Encodercut(end),10);
User.Passive.MeasuredVelocity=EncoderData(Encodercut(1):Encodercut(end),5);
User.Passive.Torque=EncoderData(Encodercut(1):Encodercut(end),7);

% 2.1) AAN 1, k=10, d=2
EMGcut=find(EMGCommand==102);
User.AAN.Biceps=EMGData(EMGcut(1):EMGcut(end),1);
User.AAN.Triceps=EMGData(EMGcut(1):EMGcut(end),2);

Encodercut=find(EncoderData(:,1)==102);
User.AAN.MeasuredPosition=EncoderData(Encodercut(1):Encodercut(end),4);
User.AAN.DesiredPosition=EncoderData(Encodercut(1):Encodercut(end),10);
User.AAN.MeasuredVelocity=EncoderData(Encodercut(1):Encodercut(end),5);
User.AAN.Torque=EncoderData(Encodercut(1):Encodercut(end),7);

% 2.2) AAN 2, k=5, d=1
EMGcut=find(EMGCommand==11);
User.AAN2.Biceps=EMGData(EMGcut(1):EMGcut(end),1);
User.AAN2.Triceps=EMGData(EMGcut(1):EMGcut(end),2);

Encodercut=find(EncoderData(:,1)==11);
User.AAN2.MeasuredPosition=EncoderData(Encodercut(1):Encodercut(end),4);
User.AAN2.DesiredPosition=EncoderData(Encodercut(1):Encodercut(end),10);
User.AAN2.MeasuredVelocity=EncoderData(Encodercut(1):Encodercut(end),5);
User.AAN2.Torque=EncoderData(Encodercut(1):Encodercut(end),7);

% 3) Antig
EMGcut=find(EMGCommand==103);
User.Antig.Biceps=EMGData(EMGcut(1):EMGcut(end),1);
User.Antig.Triceps=EMGData(EMGcut(1):EMGcut(end),2);

Encodercut=find(EncoderData(:,1)==103);
User.Antig.MeasuredPosition=EncoderData(Encodercut(1):Encodercut(end),4);
User.Antig.DesiredPosition=EncoderData(Encodercut(1):Encodercut(end),10);
User.Antig.MeasuredVelocity=EncoderData(Encodercut(1):Encodercut(end),5);
User.Antig.Torque=EncoderData(Encodercut(1):Encodercut(end),7);

% 4) Transparent
EMGcut=find(EMGCommand==104);
User.Transparent.Biceps=EMGData(EMGcut(1):EMGcut(end),1);
User.Transparent.Triceps=EMGData(EMGcut(1):EMGcut(end),2);

Encodercut=find(EncoderData(:,1)==104);
User.Transparent.MeasuredPosition=EncoderData(Encodercut(1):Encodercut(end),4);
User.Transparent.DesiredPosition=EncoderData(Encodercut(1):Encodercut(end),10);
User.Transparent.MeasuredVelocity=EncoderData(Encodercut(1):Encodercut(end),5);
User.Transparent.Torque=EncoderData(Encodercut(1):Encodercut(end),7);

% 5) Resistive
EMGcut=find(EMGCommand==105);
User.Resistive.Biceps=EMGData(EMGcut(1):EMGcut(end),1);
User.Resistive.Triceps=EMGData(EMGcut(1):EMGcut(end),2);

Encodercut=find(EncoderData(:,1)==105);
User.Resistive.MeasuredPosition=EncoderData(Encodercut(1):Encodercut(end),4);
User.Resistive.DesiredPosition=EncoderData(Encodercut(1):Encodercut(end),10);
User.Resistive.MeasuredVelocity=EncoderData(Encodercut(1):Encodercut(end),5);
User.Resistive.Torque=EncoderData(Encodercut(1):Encodercut(end),7);

% 6) Challenging
EMGcut=find(EMGCommand==106);
User.Challenging.Biceps=EMGData(EMGcut(1):EMGcut(end),1);
User.Challenging.Triceps=EMGData(EMGcut(1):EMGcut(end),2);

% User.Challenging.Biceps=EMGData(701513:EMGcut(end),1);
% User.Challenging.Triceps=EMGData(701513:EMGcut(end),2);

Encodercut=find(EncoderData(:,1)==106);
User.Challenging.MeasuredPosition=EncoderData(Encodercut(1):Encodercut(end),4);
User.Challenging.DesiredPosition=EncoderData(Encodercut(1):Encodercut(end),10);
User.Challenging.MeasuredVelocity=EncoderData(Encodercut(1):Encodercut(end),5);
User.Challenging.Torque=EncoderData(Encodercut(1):Encodercut(end),7);

% 7) Smartini=Scalini di marta to validate the anti-g
EMGcut=find(EMGCommand==10);
User.Smartini.Biceps=EMGData(EMGcut(1):EMGcut(end),1);
User.Smartini.Triceps=EMGData(EMGcut(1):EMGcut(end),2);

Encodercut=find(EncoderData(:,1)==10);
User.Smartini.MeasuredPosition=EncoderData(Encodercut(1):Encodercut(end),4);
User.Smartini.DesiredPosition=EncoderData(Encodercut(1):Encodercut(end),10);
User.Smartini.MeasuredVelocity=EncoderData(Encodercut(1):Encodercut(end),5);
User.Smartini.Torque=EncoderData(Encodercut(1):Encodercut(end),7);

% Name (to be commented if needed)
User.Name = filename(1:8);

figure,
plot(User.AAN2.Biceps)
hold on,plot(User.AAN2.MeasuredPosition)

% save('User3-Marta.mat','User')


% User.Smartini.MeasuredVelocity(1:58756)=[];