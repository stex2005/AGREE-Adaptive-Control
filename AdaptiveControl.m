close all
clear all
clc

EncoderData=readtable('AGREE-Manager-2020-08-05-16-36-14.csv');


% Column 1: elapsed time
% Column 3: position measured
% Column 4: velocity measured
% Column 7: J status
% Column 16: position desired
% Column 17: Stiffness
% Column 18: damping


EncoderData.Var1=EncoderData.Var1./1000;
% elapsed_time_s=AGREE(:,1)./1000;
% EncoderData.VelocityDesired=diff(EncoderData.Var16).*10;%%
EncoderData.Var17=EncoderData.Var17./1000;
EncoderData.Var18=EncoderData.Var18./1000;
duration=4; %AGREE(:,9);
window_length=duration*10;


%% Filtering of position and velocity desired
fc = 7;
fs = 100;

[b,a] = butter(1,fc/(fs/2));

EncoderData.PositionDesiredFilt = filtfilt(b,a,EncoderData.Var16);
VelocityDesired=diff(EncoderData.PositionDesiredFilt).*10;

% time=[1:0.01:40];
% fifthorder=10*(time/40).^3-15*(time/40).^4+6*(time/40).^5;
% velocityfifth=diff(fifthorder);

figure, plot(EncoderData.Var3)
hold on
%plot(EncoderData.PositionDesiredFilt)
%hold on
plot(EncoderData.Var7)
hold on
plot(EncoderData.Var17)
hold on
plot(EncoderData.Var18)
legend('Position measured','Position desired','Status','Stiffness','Damping')


%% Correct for the delay of the measured position
% correction_factor=1;
% position_measured(correction_factor,:)=[];
% position_desired([size(position_measured,1)+correction_factor:end],:)=[];
% elapsed_time_s([size(position_measured,1)+correction_factor:end],:)=[];
% VelocityDesired([size(position_measured,1)+correction_factor:end],:)=[];
% velocity_measured([size(position_measured,1)+correction_factor:end],:)=[];
% stiffness_normalized([size(position_measured,1)+correction_factor:end],:)=[];
% damping_normalized([size(position_measured,1)+correction_factor:end],:)=[];
% status([size(position_measured,1)+correction_factor:end],:)=[];

% %% Plot Impedance Control Status
% 
% figure();
% i=0;
% 
% i_max = 2;
% i=i+1;
% subplot(i_max,1,i)
% plot(elapsed_time_s(1:end-1),position_measured(1:end-1))
% hold on
% plot(elapsed_time_s(1:end-1),position_desired(1:end-1))
% %hold on
% %plot(elapsed_time_s(1:end-1),position_desired_filt)
% 
% legend('Measured Position [rad]','Desired Position [rad]');
% 
% % i=i+1;
% % subplot(i_max,1,i)
% % plot(elapsed_time_s(1:end-2),velocity_measured(1:end-2));
% % %hold on
% % %plot(elapsed_time_s(1:end-2),velocity_desired(1:end-1));
% %  hold on
% %  plot(elapsed_time_s(1:end-2),velocity_desired_filt);
% %  legend('Measured Velocity [rad]','Desired Velocity [rad]','Velocity filt');
% % 
% i=i+1;
% subplot(i_max,1,i)
% plot(elapsed_time_s,stiffness)
% hold on
% plot(elapsed_time_s,damping)
% %hold on
% %plot(elapsed_time_s,status)
% legend('Stiffness','Damping');

%% Starting point of setpoint

StartingPoint=find(EncoderData.Var7==7);
StartingPoint=StartingPoint(1)+20;

% for i=3:length(VelocityDesired)
%     if (VelocityDesired(i)~=0 && EncoderData.PositionDesiredFilt(i)>0 && EncoderData.PositionDesiredFilt(i)~=EncoderData.PositionDesiredFilt(i-2))
%    starting_point=i;
%         break;
%     end
% end

%% Cut the initial part without setpoint
EncoderData(1:StartingPoint,:)=[];
VelocityDesired(1:StartingPoint,:)=[];


    
%% Find the kinematic onset of movement
%[Onset_time,Onset_discrete]=FindOnset(EncoderData.PositionDesiredFilt,1);
Onset=15;

%% Plot Impedance Control Find

% figure();
% i=0;
% 
% 
% i_max = 2;
% i=i+1;
% subplot(i_max,1,i)
% plot(elapsed_time_s,position_measured)
% hold on
% plot(elapsed_time_s,position_desired)
% 
% legend('Measured Position [rad]','Desired Position [rad]');
% xlabel('Time [Sec]')
% ylabel('Position [Rad]')
% set(gca, 'FontName', 'Times New Roman','Fontsize',14)

% i=i+1;
% subplot(i_max,1,i)
% plot(elapsed_time_s,velocity_measured);
% hold on
% plot(elapsed_time_s,velocity_desire));
% % plot(elapsed_time_s,AGREE.I_position_des_rad-AGREE.J_position_rad
% hold on
% legend('Measured Velocity [rad/sec]','Desired Velocity [rad/sec]');


% i=i+1;
% subplot(i_max,1,i)
% yyaxis left
% plot(elapsed_time_s,stiffness_normalized)
% ylabel ('Stiffness')
% yyaxis right
% hold on
% plot(elapsed_time_s,damping_normalized)
% ylabel ('Damping')
% set(gca, 'FontName', 'Times New Roman')
% set(gca, 'FontName', 'Times New Roman','Fontsize',14)
% 
% %plot(elapsed_time_s,status)
% legend('Stiffness','Damping');
% xlabel('Time [Sec]')


%% Adaptive Control
StiffnessNormalized=EncoderData.Var17/max(EncoderData.Var17);

% Adaptive law: alpha(k+1)= f*alpha(k) + f*P(k)
% 0<alpha<1
% f+g=1, f=forgetting factor, g=gain factor
f=0.9;
g=0.1;
% P=performance measurement, 0<=P<=1
lambda=0.01;
PositionError=abs(EncoderData.PositionDesiredFilt-EncoderData.Var3);
PositionError_filtered=0;
PositionError_filtered=lambda*(PositionError)+(1-lambda)*PositionError_filtered;


VelocityError=VelocityDesired-EncoderData.Var4(1:end-1,:);
VelocityError_filtered=0;
VelocityError_filtered=lambda*(VelocityError)+(1-lambda)*VelocityError_filtered;

weight1=1;
weight2=0;
P=weight1*PositionError_filtered(1:end-1,:)+weight2*VelocityError_filtered;

% P=P/max(P);

alpha=ones(size(EncoderData,1)-1,1);
ErrorThreshold=0.02;
AdaptiveStiffness=EncoderData.Var17;

for i=2:size(EncoderData,1)-1
    if (P(i)>ErrorThreshold)
        alpha(i)=f.*alpha(i-1)+g.*P(i);
        AdaptiveStiffness(i)=round(alpha(i).*EncoderData.Var17(i-1));
    else
        AdaptiveStiffness(i)=round(EncoderData.Var17(i-1));
    end
end

% for i=2:size(EncoderData,1)-1
%     if (P(i)>ErrorThreshold)
%         AdaptiveStiffness(i)=round(AdaptiveStiffness(i-1)+g.*P(i));
%     else
%         AdaptiveStiffness(i)=round(EncoderData.Var17(i-1));
%     end
% end

figure, yyaxis left, plot(PositionError_filtered)
hold on, yyaxis right, plot(EncoderData.Var17)
hold on, plot(AdaptiveStiffness,'g')

% Metodo 2 - Pellivan
alpha=(PositionError_filtered-min(PositionError_filtered))/(max(PositionError_filtered)-min(PositionError));
kd_min=5;
kd_max=30;
A=round((1-alpha)*kd_min+alpha*kd_max);
tau=40;
AdaptiveStiffness=EncoderData.Var17;
for i=2:size(EncoderData,1)-1
    AdaptiveStiffness(i)=((1-1/tau)*AdaptiveStiffness(i-1)+A(i)/tau);
end
AdaptiveStiffness=round(AdaptiveStiffness);


% Metodo 3
beta=g;
gamma=f;

for i=2:size(EncoderData,1)-1
    AdaptiveStiffness(i)=AdaptiveStiffness(i-1)+beta*PositionError_filtered(i)-gamma;
end
AdaptiveStiffness=round(AdaptiveStiffness);


figure,
subplot(5,1,1),plot(EncoderData.Var3)
hold on, plot(EncoderData.PositionDesiredFilt)
legend('Position measured','Position desired')
title('Position')

subplot(5,1,2),plot(PositionError)
hold on, plot(PositionError_filtered)
legend('Position error','Position error filtered')
title('Position error')

subplot(5,1,3),plot(VelocityError)
hold on, plot(VelocityError_filtered)
legend('Velocity error','Velocity error filtered')
title('Velocity error')

subplot(5,1,4),
yyaxis left 
plot(P)
hold on
yyaxis right
plot(EncoderData.Var17)
legend('Performance measure','Stiffness')
title('Performance')

subplot(5,1,5),
yyaxis left
plot(StiffnessNormalized)
hold on, 
yyaxis right
plot(AdaptiveStiffness)
legend('Stiffness','Adaptive stiffness')
title('Adaptation')

%% outcome measure
deviation_y=0;
counter_initiated_movement=0;
time_to_initiate_mov=window_length;
success_rate_number=0;


% Column 1: elapsed time
% Column 3: position measured
% Column 4: velocity measured
% Column 7: J status
% Column 16: position desired
% Column 17: Stiffness
% Column 18: damping

for i_sample=window_length+1:length(EncoderData.PositionDesiredFilt) %i_sample=primo campione della finestra successiva, i.e. parto da 41 se finestra è 40
i_sample
    ith_window=i_sample-window_length; % campione di partenza della i-esima finestra lunga window_length sulla quale calcolo tutte le outcome 
    
    position_window_measured=EncoderData.Var3(i_sample-window_length:i_sample-1);
    position_window_desired=EncoderData.PositionDesiredFilt(i_sample-window_length:i_sample-1);

    velocity_window_measured=EncoderData.Var4(i_sample-window_length:i_sample-1);
    velocity_window_desired=VelocityDesired(i_sample-window_length:i_sample-1);

    jerk_window_measured=diff(diff(diff(position_window_measured)));
    


% 2) LOSS OF FRACTIONATED MOVEMENT

% TASK EFFICIENCY 
    % Integral error along trajectory //
    trajectory_integral_desired=cumtrapz(position_window_desired);
    trajectory_integral_measured=cumtrapz(position_window_measured);
    integral_desired(ith_window)=trajectory_integral_desired(end);%integral under the desired trajectory (cumulative measure)
    integral_measured(ith_window)=trajectory_integral_measured(end); %integral under the measured trajectory (cumulative measure)
    absolute_error(ith_window)=abs(trajectory_integral_desired(end)-trajectory_integral_measured(end)); %Integral error between desired and measured trajectory within one window

% 3) ABNORMAL MUSCLE TONE

% RANGE OF MOTION
    % Maximum and minimum reached amplitude - Rad //
    max_reached_amplitude(ith_window)=max(position_window_measured);
    min_reached_amplitude(ith_window)=min(position_window_measured);
    
% EASE OF MOVEMENTS
 
    % Mean velocity //
    velocity_measured_window_mean(ith_window)=mean(abs(velocity_window_measured));
    velocity_desired_window_mean(ith_window)=mean(abs(velocity_window_desired));
   
    % Difference between desired and measured mean velocity
    velocity_difference(ith_window)=velocity_measured_window_mean(ith_window)-velocity_desired_window_mean(ith_window);
    
    % Peak velocity //
    velocity_window_peak(ith_window)=max(abs(velocity_window_measured));
    
    % Peak Speed Ratio //
    peak_speed_ratio(ith_window)=velocity_measured_window_mean(ith_window)/velocity_window_peak(ith_window);

    
% Useful parameters computation through for cycle - sbatta....
    initiated_movement=0;
    time_to_initiate_mov=0;
    time_to_terminate_mov=0;
    time_to_initiate_mov_found=0;
    time_to_terminate_mov_found=1;
    counter_mapr=0;
    
        for i=1:window_length-1
            
    % Ability to initiate movement //
            if (i < 0.3*window_length && initiated_movement==0 && abs(velocity_window_measured(i)) > 0.2*velocity_window_peak(ith_window)) 
                initiated_movement=1;
               % counter_initiated_movement=counter_initiated_movement+1 %Con il banco di prova dove non finestro in task fa schifo
            end
            
    % Movement planning - Time to peak speed //
            if (abs(velocity_window_measured(i)) == velocity_window_peak(ith_window))
                time_to_peakspeed(ith_window)=i;
            end
    % Temporal efficiency - Time between movement onset and termination //
            if (time_to_initiate_mov_found==0 && abs(velocity_window_measured(i)) > 0.2*velocity_window_peak(ith_window))
                time_to_initiate_mov=i;
                time_to_initiate_mov_found=1;
                time_to_terminate_mov_found=0;
            end       
            if (time_to_terminate_mov_found==0 && i > time_to_initiate_mov && abs(velocity_window_measured(i)) < 0.2*velocity_window_peak(ith_window))
                time_to_terminate_mov=i;
                time_to_terminate_mov_found=1;
            end 
           
            if (abs(velocity_window_measured(i)) > (0.2*velocity_window_peak))
                counter_mapr=counter_mapr+1;
            end
           
            
            frac_length_dp_measured(i)=sqrt((position_window_measured(i+1)-position_window_measured(i))^2+(1)^2);
            frac_length_dp_desired(i)=sqrt((position_window_desired(i+1)-position_window_desired(i))^2+(1)^2);
        end
    
    
% SMOOTHNESS
    % Speed difference //
    speed_difference(ith_window)=max(abs(velocity_window_measured))-min(abs(velocity_window_measured));

    % Mean acceleration 
    acceleration_measured_mean(ith_window)=mean(diff(velocity_window_measured));
    acceleration_desired_mean(ith_window)=mean(diff(velocity_window_desired));

    % Peak acceleration ratio 
    acceleration_ratio(ith_window)=acceleration_measured_mean(ith_window)/max(diff(velocity_window_measured));
    
    % Mean Arrest Period Ratio //
    mapr(ith_window)=counter_mapr/window_length;
    
    % Normalized Jerk (dimensioneless) - Germanotta et al 2015 
    normalized_jerk(ith_window)=sqrt(0.5*(window_length^5/(sum(frac_length_dp_measured))^2)*(integral_measured(ith_window))^2);

    % Dimensional Jerk - Rohrer et al 2002 
    jerk_measured_mean=mean(jerk_window_measured);
    not_normalized_jerk(ith_window)=-jerk_measured_mean/velocity_window_peak(ith_window);
    
    % Spectral arc length - Balasubramanianet al 2012
    %Ts=0.01; %Sampling time in seconds
    %spectral_arc_length(ith_window)=SpectralArcLength(velocity_window_measured,Ts);

    
% 4) LOSS OF SOMATOSENSATION

% MOVEMENT PLANNING
%     % Percentage of initiated movements - Mazzoleni 2011 -> NON VA BENE
%     CON CICLI CONTINUI COME BANCO DI PROVA
%     percentage_initiated_movement=counter_initiated_movement*100/(length(position_measured));
    
    % Aiming angle - Rad //
    starting_point=[1; position_window_desired(1)];
    ending_point_desired=[window_length; position_window_desired(window_length)];
    ending_point_measured=[window_length; position_window_measured(window_length)];

    vector_desired=ending_point_desired-starting_point
    vector_measured=ending_point_measured-starting_point
    aiming_angle(ith_window)=acosd(dot(vector_desired,vector_measured)/(norm(vector_desired)*norm(vector_measured)));
   
    
% TEMPORAL EFFICIENCY
    % Time between movement onset and termination //
    movement_duration(ith_window)=time_to_terminate_mov-time_to_initiate_mov;

% ACCURACY
    % Mean deviation from desired trajectory //
    deviation_y=position_window_measured-position_window_desired;
    deviation_y_mean(ith_window)=sum(deviation_y)/window_length;
    
    % Normalized path length //
    trajectory_measured_length(ith_window)=sum(frac_length_dp_measured);
    trajectory_desired_length(ith_window)=sum(frac_length_dp_desired);
    normalized_path_length(ith_window)=trajectory_measured_length/trajectory_desired_length;

    % End point error - Schwarz et al 2020 //
    end_point_error_ith(ith_window)=ending_point_desired(2)-ending_point_measured(2);
    end_point_error=mean(end_point_error_ith);

% PRECISION
    % Variable error - Schwarz et al 2020
    variable_error=std(end_point_error_ith);
    
% EFFICACY
    % Success rate - Schwarz et al 2020 //
    if (abs(end_point_error_ith(ith_window)) < 0.10*trajectory_desired_length)
        success_rate_number=success_rate_number+1;
    end
    
    success_rate=success_rate_number*100/(i_sample-window_length);
    
    
    clear position_window_measured;
    clear position_window_desired;
    clear deviation_y;
    clear velocity_window_measured;
    clear velocity_window_desired;
    clear jerk_desired;
    clear jerk_measured;
    clear length_dp_measured;
    clear length_dp_desired;
end





%% Position dependent measures

figure();
i=0;

i_max = 2;
i=i+1;
subplot(i_max,1,i)
plot(elapsed_time_s(1:end-window_length),position_measured(1:end-window_length))
hold on
plot(elapsed_time_s(1:end-window_length),position_desired(1:end-window_length))
legend('Measured Position [rad]','Desired Position [rad]');

i=i+1;
subplot(i_max,1,i)
plot(elapsed_time_s(1:end-window_length),integral_measured);
hold on
plot(elapsed_time_s(1:end-window_length),integral_desired);
hold on
plot(elapsed_time_s(1:end-window_length),absolute_error);
legend('Integral error measured','Integral error desired','Absolute error');


% Mean deviation - Integral error
figure();
i=0;

i_max = 3;
i=i+1;
subplot(i_max,1,i)
plot(elapsed_time_s(1:end-window_length),position_desired(1:end-window_length))
hold on
plot(elapsed_time_s(1:end-window_length),position_measured(1:end-window_length))
hold on
legend ('Measured Position [rad]','Desired Position [rad]')

i=i+1;
subplot(i_max,1,i)
plot(elapsed_time_s(1:end-window_length),deviation_y_mean)
hold on
plot(elapsed_time_s(1:end-window_length),absolute_error/10)
hold on
xL = xlim;
line(xL, [0 0],'Color','k','LineStyle','--','LineWidth',1);  %y-axis
legend ('Mean Deviation on Y axes','Integral error along trajectory','Zero')

i=i+1;
subplot(i_max,1,i)
plot(elapsed_time_s,stiffness_normalized)
hold on
plot(elapsed_time_s,damping_normalized)
plot(elapsed_time_s,status)
legend('Stiffness','Damping','Mode');



% i=i+1;
% subplot(i_max,1,i)
% plot(elapsed_time_s(1:end-window_length),max_reached_amplitude)
% hold on
% plot(elapsed_time_s(1:end-window_length),min_reached_amplitude)
% hold on
% plot(elapsed_time_s(1:end-window_length),position_measured(1:end-window_length))
% legend ('Maximum amplitude [Rad]','Minimum amplitude [Rad]','Position measured [Rad]')


% Aiming angle
figure();
i=0;

i_max = 3;

i=i+1;
subplot(i_max,1,i)
plot(elapsed_time_s(1:end-window_length),position_desired(1:end-window_length))
hold on
plot(elapsed_time_s(1:end-window_length),position_measured(1:end-window_length))
hold on
plot(elapsed_time_s(1:end-window_length),aiming_angle)
legend ('Desired Position [rad]','Measured Position [rad]','Aiming angle [Rad]')


i=i+1;
subplot(i_max,1,i)
plot(elapsed_time_s(1:end-window_length),deviation_y_mean)
hold on
plot(elapsed_time_s(1:end-window_length),aiming_angle)
hold on
plot(elapsed_time_s(1:end-window_length),absolute_error/10)
legend ('Mean Deviation on Y axes [rad]','Aiming angle [Rad]','Integral error along trajectory')

i=i+1;
subplot(i_max,1,i)
plot(elapsed_time_s,stiffness_normalized)
hold on
plot(elapsed_time_s,damping_normalized)
plot(elapsed_time_s,status)
legend('Stiffness','Damping','Mode');



% Normalized trajectory length
figure();
i=0;

i_max = 4;

i=i+1;
subplot(i_max,1,i)
plot(elapsed_time_s(1:end-window_length),position_desired(1:end-window_length))
hold on
plot(elapsed_time_s(1:end-window_length),position_measured(1:end-window_length))
legend ('Desired Position [rad]','Measured Position [rad]')


i=i+1;
subplot(i_max,1,i)
plot(elapsed_time_s(1:end-window_length),trajectory_desired_length)
hold on
plot(elapsed_time_s(1:end-window_length),trajectory_measured_length)
legend ('Trajectory desired length','Trajectory measured length',' Normalized path length')

i=i+1;
subplot(i_max,1,i)
plot(elapsed_time_s(1:end-window_length),normalized_path_length)
legend('Normalized path length')

i=i+1;
subplot(i_max,1,i)
plot(elapsed_time_s,stiffness_normalized)
hold on
plot(elapsed_time_s,damping_normalized)
plot(elapsed_time_s,status)
legend('Stiffness','Damping','Mode');



i=i+1;
subplot(i_max,1,i)
plot(elapsed_time_s(1:end-window_length),peak_speed_ratio)
legend ('Peak Speed Ratio')




i=i+1;
subplot(i_max,1,i)
plot(elapsed_time_s(1:end-window_length),normalized_path_length)
legend ('Normalized path length')


%% Time dependent measures

figure();
i=0;

i_max = 4;

i=i+1;
subplot(i_max,1,i)
plot(elapsed_time_s(1:end-window_length),position_desired(1:end-window_length))
hold on
plot(elapsed_time_s(1:end-window_length),position_measured(1:end-window_length))
legend ('Measured Position [rad]','Desired Position [rad]')


i=i+1;
subplot(i_max,1,i)
plot(elapsed_time_s(1:end-window_length),time_to_peakspeed)
legend ('Time to peak speed')

i=i+1;
subplot(i_max,1,i)
plot(elapsed_time_s(1:end-window_length),movement_duration)
legend('Movement duration')

i=i+1;
subplot(i_max,1,i)
plot(elapsed_time_s,stiffness_normalized)
hold on
plot(elapsed_time_s,damping_normalized)
plot(elapsed_time_s,status)
legend('Stiffness','Damping','Mode');


%% Velocity dependent measures

figure();
i=0;

i_max = 4;

i=i+1;
subplot(i_max,1,i)
plot(elapsed_time_s(1:end-window_length),position_desired(1:end-window_length))
hold on
plot(elapsed_time_s(1:end-window_length),position_measured(1:end-window_length))
legend ('Measured Position [rad]','Desired Position [rad]')


i=i+1;
subplot(i_max,1,i)
plot(elapsed_time_s(1:end-window_length),velocity_measured_window_mean)
hold on
plot(elapsed_time_s(1:end-window_length),velocity_measured_window_mean-velocity_desired_window_mean)
hold on
plot(elapsed_time_s(1:end-window_length),peak_speed_ratio)
legend ('Mean Velocity','Mean velocity difference','Peak speed ratio')

i=i+1;
subplot(i_max,1,i)
plot(elapsed_time_s(1:end-window_length),mapr)
legend('Mean Arrest Period Ration')

i=i+1;
subplot(i_max,1,i)
plot(elapsed_time_s,stiffness_normalized)
hold on
plot(elapsed_time_s,damping_normalized)
plot(elapsed_time_s,status)
legend('Stiffness','Damping','Mode');


%% Acceleration dependent measures

figure();
i=0;

i_max = 4;

i=i+1;
subplot(i_max,1,i)
plot(elapsed_time_s(1:end-window_length),position_desired(1:end-window_length))
hold on
plot(elapsed_time_s(1:end-window_length),position_measured(1:end-window_length))
legend ('Measured Position [rad]','Desired Position [rad]')


i=i+1;
subplot(i_max,1,i)
plot(elapsed_time_s(1:end-window_length),acceleration_measured_mean-acceleration_desired_mean)
legend ('Mean acceleration difference')

i=i+1;
subplot(i_max,1,i)
plot(elapsed_time_s(1:end-window_length),acceleration_ratio)
legend ('Peak acceleration ratio')


i=i+1;
subplot(i_max,1,i)
plot(elapsed_time_s,stiffness_normalized)
hold on
plot(elapsed_time_s,damping_normalized)
plot(elapsed_time_s,status)
legend('Stiffness','Damping','Mode');

%% Jerk

figure();
i=0;

i_max = 4;

i=i+1;
subplot(i_max,1,i)
plot(elapsed_time_s(1:end-window_length),position_desired(1:end-window_length))
hold on
plot(elapsed_time_s(1:end-window_length),position_measured(1:end-window_length))
legend ('Measured Position [rad]','Desired Position [rad]')


i=i+1;
subplot(i_max,1,i)
plot(elapsed_time_s(1:end-window_length),normalized_jerk)
legend ('Normalized Jerk')

i=i+1;
subplot(i_max,1,i)
plot(elapsed_time_s(1:end-window_length),not_normalized_jerk)
legend ('Not normalized Jerk')


i=i+1;
subplot(i_max,1,i)
plot(elapsed_time_s,stiffness_normalized)
hold on
plot(elapsed_time_s,damping_normalized)
plot(elapsed_time_s,status)
legend('Stiffness','Damping','Mode');