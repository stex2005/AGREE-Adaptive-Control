function [emg_processed] = EMGFilter(emg_raw,ft_l,ft_h,fc)

% 0) 3th order Butterworth filter

% High Pass Filter coefficient definition
Wn_h= ft_h/(fc/2);
[b_h,a_h] = butter(3,Wn_h,'high');

% Low Pass Filter coefficient definition
Wn_l= ft_l/(fc/2);
[b_l,a_l] = butter(3,Wn_l,'low'); 

% 1) High Pass filter
emg_high=filtfilt(b_h,a_h,emg_raw);

% 2) Rectification
emg_rect=abs(emg_high);

% 3) Low Pass filter
emg_low=filtfilt(b_l,a_l,emg_rect);

emg_processed=emg_low;
% 4) Normalization
% emg_processed=(emg_low)/(Maximum);



end

