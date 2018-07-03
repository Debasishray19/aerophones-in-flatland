function audioGenfunc (audioData)

    % Define audio sampling rate, Frame width & buffer size
    Fs = 44100; % Number of samples per second
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
        nfft2 = 2^nextpow2(nfft);

        fft_signal = fft(audioData, nfft2);
        fft_signal = fft_signal(1:nfft2/2);
        fft_signal = fft_signal/max(fft_signal);
        xfft = Fs*(0:nfft2/2 -1)/nfft2;
        figure;
        %subplot(2,1,1); plot(dt,audioData, '-r'); ylim([-3 3]);
        subplot(2,1,2); plot(xfft, abs(fft_signal), '-b'); xlim([0 44100]);

    end
    % Audio streaming loop
    rangeDefiner=0;

    for loopCount = 1:samplesPerFrame-1
        xFrameWindow = ((1+(rangeDefiner*FrameWidth)):(rangeDefiner+1)*FrameWidth);
        signal = audioData(xFrameWindow);
        deviceWriter(signal);
        rangeDefiner= rangeDefiner+1;
    end   
end