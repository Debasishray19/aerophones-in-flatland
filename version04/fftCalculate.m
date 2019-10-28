function [signal_fft_final] = fftCalculate(signal)
    % fftCalculate dunction computes and returns the fft of a signal

    % Signal length
    signalLen = length(signal);
    
    % Compute fft
    signal_fft = fft(signal./max(signal));
    signal_fft_abs = abs(signal_fft/signalLen);
    
    if mod(signalLen,2)==0
        signal_fft_final = signal_fft_abs(1:signalLen/2);
    else   
        signal_fft_final = signal_fft_abs(1:(signalLen-1)/2);
    end
end