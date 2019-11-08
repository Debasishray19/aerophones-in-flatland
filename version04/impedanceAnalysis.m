close all; 
clear;
rho = 1.140;
c   = 350;
fs = 15*44100;
load('impedanceData_20_ms.mat');
% load('impedanceData.mat');
% load('impedanceData_50_ms.mat');
% load('impedanceData_1_2_sec.mat');
Pr_fft = fft(Pr_Audio);
Vx_fft = fft(Vx_Vel);

Vx_fft_real = real(Vx_fft);
Vx_fft_imag = imag(Vx_fft);

Vx_theoritical_real = real(Pr_fft)./(rho*c);
Vx_theoritical_imag = imag(Pr_fft)./(rho*c);
Vx_theoritical_fft  = Vx_theoritical_real + (1i.*Vx_theoritical_imag);

% Calculating Z impedance
z_actual = Pr_fft./Vx_fft;
z_theoritical = Pr_fft./Vx_theoritical_fft;

% Source Signal samples
% plotFrequencyPhase(excitationV(1:length(Pr_Audio)), fs)
    
% Plot resistance
figure;  hold on;
plot(real(z_actual));
plot(real(z_theoritical));
xlabel('frequency in Hz');
ylabel('Resistance');
axis 'auto y';
yLim = ylim();
axis([2 23050 yLim(1) yLim(2)]);
legend('zActual','zTheoritical');
hold off;

% Plot reactance
figure; hold on;
plot(imag(z_actual)); 
plot(imag(z_theoritical));
xlabel('frequency in Hz');
ylabel('Reactance');
axis 'auto y';
yLim = ylim();
axis([2 23050 yLim(1) yLim(2)])
legend('zActual','zTheoritical');
hold off;

% Plot fft(Vx_Vel) real
figure; hold on;
plot(real(Vx_fft)); 
plot(real(Vx_theoritical_fft));
xlabel('frequency in Hz');
ylabel('Velocity in m/s');
axis 'auto y';
yLim = ylim();
axis([2 23050 yLim(1) yLim(2)])
legend('VxActual','VxTheoritical');
hold off;

% Plot fft(Vx_Vel) imag
figure;hold on;
plot(imag(Vx_fft)); 
plot(imag(Vx_theoritical_fft));
xlabel('frequency in Hz');
ylabel('Velocity in m/s');
axis 'auto y';
yLim = ylim();
axis([2 23050 yLim(1) yLim(2)])
legend('VxActual','VxTheoritical');
hold off;
