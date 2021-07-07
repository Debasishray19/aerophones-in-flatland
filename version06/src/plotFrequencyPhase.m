function plotFrequencyPhase(audioSignal, sourceSignal, srate)
    
    % Define sampling rate
    fs = srate;
    audiosignalfft = fft(audioSignal);
    sourceSignalfft = fft(sourceSignal);
    
    TransferFunction = audiosignalfft./sourceSignalfft;
    finalTransferFunction = 10*log10(abs(TransferFunction).^2);
        
    audiosignalfft_phase = angle(audiosignalfft);
    n = length(audioSignal)-1;  
    df = fs/n;
    f=0:df:fs;
    
    figure
    plot(f,finalTransferFunction); % plot Fourier Transform
    title('Amplitude Spectrum Analysis');
    xlabel('Frequency [Hz]');
    ylabel('Magnitude [dB]');
    axis 'auto y';
    yLim = ylim();
    axis([2 10050 yLim(1) yLim(2)])
    
    figure
    plot(f,audiosignalfft_phase); % plot Fourier Transform
    title('Phase Spectrum Analysis');
    xlabel('Frequency [Hz]');
    ylabel('Phase [Radian]');
    axis 'auto y';
    yLim = ylim();
    axis([2 22050 yLim(1) yLim(2)])
end