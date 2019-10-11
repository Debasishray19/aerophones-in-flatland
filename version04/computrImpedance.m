function computrImpedance(pressure, velocity)

    signalLen = length(pressure);   % Compute signal length
    nfft2 = 2^nextpow2(signalLen);  % Find the next number of signalLen which is power of 2
    
    % Compute fft transformation for both pressure and velocity signal
    pressure_fft = fft (pressure, nfft2);
    velocity_fft = fft(velocity, nfft2);
    
    % To remove the mirror outputs
    pressure_fft_Final =  pressure_fft(1:nfft2/2);
    velocity_fft_Final = velocity_fft(1:nfft2/2);
    
    % Calculate impedance
    zImpedance = pressure_fft_Final./velocity_fft_Final;
    
    % Plot the impedance
    figure;
    plot(real(zImpedance));
    xlim([0 10000]);
    
    figure;
    plot(imag(zImpedance));
    xlim([0 10000]);
end