close all
clearvars -except AGREE
clc

%% Import data

% [filename, pathname] = uigetfile('*.csv',['Select a file: ']);
% if filename == 0
%     disp('No file input. Stopping simulation.');
% else
%     AGREE = importdata(filename);
% end
% 
% for i = 1:length(AGREE.data)
%     AGREE.J_elapsed_time_ms = AGREE.data(:,1);
%     AGREE.J_position_rad    = AGREE.data(:,3);
%     AGREE.J_position_des_rad= AGREE.data(:,1);
%     AGREE.J_elapsed_time_ms = AGREE.data(:,1);
%     AGREE.J_elapsed_time_ms = AGREE.data(:,1);
% end
% 

%% Plot Impedance Control Status

figure();
i=0;

elapsed_time_s = AGREE.J_elapsed_time_ms./1000;

i_max = 3;
i=i+1;
subplot(i_max,1,i)
plot(elapsed_time_s,AGREE.J_position_rad)
hold on
plot(elapsed_time_s,AGREE.I_position_des_rad)

legend('Measured Position [rad]','Desired Position [rad]');

i=i+1;
subplot(i_max,1,i)
plot(elapsed_time_s,AGREE.J_velocity_rad_s);
hold on

fc = 10;
fs = 100;

[b,a] = butter(4,fc/(fs/2));

I_position_des_rad_filt = filtfilt(b,a,AGREE.I_position_des_rad(1:end-1));

plot(elapsed_time_s(1:end-1),hampel(diff(AGREE.I_position_des_rad).*10));
% plot(elapsed_time_s,AGREE.I_position_des_rad-AGREE.J_position_rad
hold on

legend('Measured Velocity [rad/s]','Desired Velocity [rad/s]');

i=i+1;
subplot(i_max,1,i)
plot(elapsed_time_s,AGREE.I_stiffness/1000)
hold on
plot(elapsed_time_s,AGREE.I_damping/1000)
plot(elapsed_time_s,AGREE.J_status)
legend('Stiffness','Damping','Mode');

%% Plot Velocity 

figure();
plot(elapsed_time_s,AGREE.J_velocity_rad_s);
hold on;
plot(elapsed_time_s,AGREE.J_velocity_computed*1000);
title('Joint Velocity [rad/s]');
ylabel('Angular Velocity [rad/s]');
xlabel('Time [s]');


