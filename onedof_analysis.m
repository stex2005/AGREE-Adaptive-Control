clear all
close all

%% Dependencies
%Adding functions written for this simulation
curr_path = pwd;
path = [curr_path '/functions'];
addpath(path);
path = [curr_path '/onedof_log'];
addpath(path);

%% Initialization
Ntrials   = 1;
color_des = ["b";"b"];
color_meas = ["r";"r"];
Modes = [101,102,103,104,105,106];
onedof_log = {Ntrials};

for i=1:Ntrials
    
    %% Import Data
   
    [filename, pathname] = uigetfile('*.csv',['Select a file: ']);
    if filename == 0
        disp('No file input. Stopping simulation.');
    else
        onedof_log(i) = {import_onedof_log(filename)};
    end
    
    %% Modes
    for m=1:length(Modes)
        
    start = min(find(onedof_log{i}.J_status==Modes(m)));
    stop = max(find(onedof_log{i}.J_status==Modes(m)));
    
    %% Position
    figure(1);
    subplot(2,1,1);
    hold on
%     plot(onedof_log{i}.I_0_position_des_rad(start:stop));
%     plot(onedof_log{i}.E_0_position_des_rad(start:end));
    plot(-onedof_log{i}.J_0_position_rad(start:stop));
    ylim([-pi;pi]);
    subplot(2,1,2);
    hold on
    plot(-onedof_log{i}.J_0_velocity_rad_s(start:stop));
    
    %%
    
    %% Power spectrum
%     fs = 1000;
%     data = onedof_log{i}.J_0_position_rad(start:end);
%     data_l=length(onedof_log{i}.J_0_position_rad(start:end));
%     xdft = fft(data,[],1);
%     xdft = xdft(1:data_l/2+1,:);
%     psdx = (1/(fs*data_l)) * abs(xdft).^2;
%     psdx(2:end-1,:) = 2*psdx(2:end-1,:);
% 
%     freq = 0:fs/data_l:fs/2;
% 
%     figure(2);
%     hold on
%     plot(freq,psdx(:,1));
%     title('Original Periodagram ch1');
    end
    legend('Passive','Impedance','Anti-G','Transparent','Resistive','Challenging');
end

