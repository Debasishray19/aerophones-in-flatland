function frequencyPlot(audioSignal, srate)
    
    % Define sampling rate
    fs = srate;
    audiosignalfft = 10*log10(abs(fft(audioSignal)).^2);
    n = length(audioSignal)-1; 
    
    df = fs/n;
    f=0:df:fs;

    figure
    plot(f,audiosignalfft); % plot Fourier Transform
    title('Spectrum Analysis');
    xlabel('Frequency [Hz]');
    ylabel('Magnitude [dB]');
    axis 'auto y'
    yLim = ylim();
    axis([2 22050 yLim(1) yLim(2)])
end
