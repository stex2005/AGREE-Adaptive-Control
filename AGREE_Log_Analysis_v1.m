close all
clearvars -except AGREE
clc

%% Import data

[filename, pathname] = uigetfile('*.csv',['Select a file: ']);
if filename == 0
    disp('No file input. Stopping simulation.');
else
    AGREE = importdata(filename);
end

for i = 1:length(AGREE.data)
    AGREE.J_elapsed_time_ms = AGREE.data(:,1);
    AGREE.J_position_rad    = AGREE.data(:,3);
    AGREE.J_position_des_rad= AGREE.data(:,1);
    AGREE.J_elapsed_time_ms = AGREE.data(:,1);
    AGREE.J_elapsed_time_ms = AGREE.data(:,1);
end
AGREE.J_torque = AGREE.data(:,6);

%% Plot Impedance Control Status

figure();
i=0;

AGREE.elapsed_time_s = AGREE.J_elapsed_time_ms./1000;

i_max = 3;
i=i+1;
subplot(i_max,1,i)
plot(AGREE.elapsed_time_s,AGREE.J_position_rad)
hold on
plot(AGREE.elapsed_time_s,AGREE.I_position_des_rad)

legend('Measured Position [rad]','Desired Position [rad]');

i=i+1;
subplot(i_max,1,i)
plot(AGREE.elapsed_time_s,AGREE.J_velocity_rad_s);
hold on

fc = 10;
fs = 100;

[b,a] = butter(4,fc/(fs/2));

I_position_des_rad_filt = filtfilt(b,a,AGREE.I_position_des_rad(1:end-1));

plot(AGREE.elapsed_time_s(1:end-1),hampel(diff(AGREE.I_position_des_rad).*10));
% plot(elapsed_time_s,AGREE.I_position_des_rad-AGREE.J_position_rad
hold on

legend('Measured Velocity [rad/s]','Desired Velocity [rad/s]');

i=i+1;
subplot(i_max,1,i)
plot(AGREE.elapsed_time_s,AGREE.I_stiffness/1000)
hold on
plot(AGREE.elapsed_time_s,AGREE.I_damping/1000)
plot(AGREE.elapsed_time_s,AGREE.J_status)
legend('Stiffness','Damping','Mode');

%% Plot Velocity 

figure();
plot(AGREE.elapsed_time_s,AGREE.J_velocity_rad_s);
hold on;
plot(AGREE.elapsed_time_s,AGREE.J_velocity_computed*1000);
title('Joint Velocity [rad/s]');
ylabel('Angular Velocity [rad/s]');
xlabel('Time [s]');

%% Plot Transparency

figure();
t_end = 170;
yyaxis right
plot(AGREE.elapsed_time_s(1:t_end),AGREE.J_velocity_rad_s(1:t_end),'LineWidth',1);
hold on
%     plot(AGREE.elapsed_time_s(1:i),AGREE.I_position_des_rad(1:i));
ylim([-2 +2])
% xlim([0 100])
ylabel('Velocity [rad/s]','FontSize',14);
yyaxis left

plot(AGREE.elapsed_time_s(1:t_end),AGREE.J_torque_loadcell(1:t_end)/1000+0.043,'LineWidth',1);
ylim([-5 +5])
ylabel('Torque [Nm]','FontSize',14);
xlabel('Time [s]','FontSize',14);
set(gca,'FontSize', 14);
legend('Velocity [rad/s]','Transparency Torque [Nm]');

%% Dynamic Plot Impedance

figure(5);
hold off
for i = 57:length(AGREE.elapsed_time_s)
    tic
    plot(AGREE.elapsed_time_s(1:i),-AGREE.J_position_rad(1:i),'-b','LineWidth',1);
    hold on
    plot(AGREE.elapsed_time_s(1:i),-AGREE.I_position_des_rad(1:i),'-r','LineWidth',1);
    ylim([-2 +2])
    ylabel('Position [rad]');
    
%     yyaxis left
%     plot(AGREE.elapsed_time_s(1:i),AGREE.J_torque_loadcell(1:i)/1000+0.043,'LineWidth',1);
%     ylim([-10 +10])
%     ylabel('Torque [Nm]','FontSize',14);
    xlabel('Time [s]');
    t1 = toc;
%     pause(0.1-t1)
end

%% Impedance Plot

figure(6);
set(0, 'DefaultAxesFontName', 'Times');
hold off
t0 = 420;
t1 = 700;
t2 = 767;
tend = 960;

subplot(2,1,1);
plot(AGREE.elapsed_time_s(t0:tend),-AGREE.J_position_rad(t0:tend),'-r','LineWidth',1);
hold on
plot(AGREE.elapsed_time_s(t0:tend),-AGREE.I_position_des_rad(t0:tend),'-b','LineWidth',1);
ylim([-1.5 +1.5])
ylabel('Position [rad]');
legend('Measured Position','Desired Position');
set(gca,'FontSize', 12,'FontName','Times');

subplot(2,1,2);
plot(AGREE.elapsed_time_s(t0:tend),AGREE.J_torque_loadcell(t0:tend)/1000+0.043,'-b','LineWidth',1);
ylim([-7 +7])
ylabel('Torque [Nm]');
xlabel('Time [s]');
set(gca,'FontSize', 12,'FontName','Times');

%% Transparency Plot

figure(7);
set(0, 'DefaultAxesFontName', 'Times');
% hold off
t0 = 1;
tend = 200;

subplot(2,1,1);
plot(AGREE.elapsed_time_s(t0:tend),-AGREE.J_position_rad(t0:tend),'-r','LineWidth',1);
hold on
plot(AGREE.elapsed_time_s(t0:tend),-AGREE.J_velocity_rad_s(t0:tend),'-b','LineWidth',1);
ylim([-2 +2])
ylabel('Arbitrary Units [a.u.]');
legend('Measured Position [rad]','Measured Velocity [rad/s]');
set(gca,'FontSize', 12,'FontName','Times');

subplot(2,1,2);
plot(AGREE.elapsed_time_s(t0:tend),AGREE.J_torque_loadcell(t0:tend)/1000+0.043,'-b','LineWidth',1);
ylim([-2 +2])
ylabel('Torque [Nm]');
xlabel('Time [s]');
set(gca,'FontSize', 12,'FontName','Times');


%%
yaxis([-2;+2])
hold on;
% plot(elapsed_time_s,AGREE.J_velocity_computed*1000);
title('Joint Velocity [rad/s]');
ylabel('Angular Velocity [rad/s]');
xlabel('Time [s]');
yyaxis right

plot(elapsed_time_s,AGREE.J_torque_loadcell+43);




