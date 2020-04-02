close all;
clc; clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEFINE TIME CONSTANTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SECOND = 1;    % Unit of time
SRATE = 44100*5;                 % Sampling rate of sound card
dt    = 1/SRATE;               % Time period
total_audio_time = 1*SECOND;   % Temporal length of sound generation
t = 0:dt:total_audio_time-dt;  % Temporal time steps
STEPS = length(t);

output_ug   = zeros(1, STEPS);
[airParam, vf_structuralParam, vf_flowParam, vf_matParam] = setVocalFoldParams();

for time_counter = 1:STEPS
    
    [vf_flowParam] = vf_TwoMassModel(SRATE, airParam, vf_structuralParam, vf_flowParam);
    
    output_ug(time_counter) = vf_flowParam.ug_next;   
end

% Plot the vocal fold volume velocity after normalization
norm_ug = output_ug./max(output_ug);
plot(t(1:15000), norm_ug(1:15000));
ylim([0 2]);
xlabel('time');
ylabel('volume velocity');