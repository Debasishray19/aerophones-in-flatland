function audioGenfunc (audioData, sourceData)
    % Initialize MATLAB
    close all;
    clc;
    
    % Define audio sampling rate, Frame width & buffer size
    Fs = 44100; % Number of samples per second
    Ts = 1/Fs;  % Sapling period
    timeA = (1:length(sourceData)).*Ts;
    
    FrameWidth = 1250;
    BufferSize = 1250;
    
    % Number of samples per frame
    samplesPerFrame = round(length(audioData)/FrameWidth);
    
    % Create system object to write to buffer
    deviceWriter = audioDeviceWriter('SampleRate', Fs,...
    'SupportVariableSizeInput',true, 'BufferSize',BufferSize);
    
    % Plot Spectrum
    plotSpectrum = true;

    if plotSpectrum
        nfft = length(audioData);
        src_nfft = length(sourceData);
        
        nfft2 = 2^nextpow2(nfft);
        src_nfft2 = 2^nextpow2(src_nfft);
        
        fft_signal = fft(audioData, nfft2);
        src_fft_signal = fft(sourceData, src_nfft2);
        
        fft_signal = fft_signal(1:nfft2/2);
        src_fft_signal = src_fft_signal(1:nfft2/2);
        
        fft_signal = fft_signal/max(fft_signal);
        src_fft_signal = src_fft_signal/max(src_fft_signal);
        
        xfft = Fs*(0:nfft2/2 -1)/nfft2;
        src_xfft = Fs*(0:src_nfft2/2 -1)/src_nfft2;
        
        figure;
        subplot(4,1,1); plot(timeA,sourceData, '-r');
        title('Source Signal');
        xlabel('Time');
        ylabel('Amplitude');
  
        subplot(4,1,2); plot(timeA,audioData, '-b');
        title('Output Signal');
        xlabel('Time');
        ylabel('Amplitude');
        
        subplot(4,1,3); plot(src_xfft, abs(src_fft_signal), '-r'); xlim([0 44100]);
        title('Source Spectrum');
        xlabel('Freequency in Hz');
        ylabel('Amplitude');
        
        subplot(4,1,4); plot(xfft, abs(fft_signal), '-b'); xlim([0 44100]);
        title('Output Spectrum');
        xlabel('Freequency in Hz');
        ylabel('Amplitude');
    end
    
    % Audio streaming loop
    rangeDefiner=0;

    for loopCount = 1:samplesPerFrame-1
        xFrameWindow = ((1+(rangeDefiner*FrameWidth)):(rangeDefiner+1)*FrameWidth);
        signal = audioData(xFrameWindow);
        deviceWriter(signal');
        rangeDefiner= rangeDefiner+1;
    end   
end