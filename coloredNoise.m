function [signal,time] = coloredNoise(pv)
% Generate nrTimePoints of colored noise, sampled at a specified sampling rate 
% and with a specified drop off of the power spectrum with frequency.
% 
% Pink Noise (power ~1/f):
% pn = coloredNoise(nrSamples = 10000,samplingRate =1e3, alpha=1)
arguments
    pv.nrTimePoints (1,1) double {mustBeInteger, mustBeNonnegative} = 1000
    pv.samplingRate (1,1) double {mustBeNonnegative} = 1
    pv.alpha (1,1) double {mustBeInteger, mustBeNonnegative} =1 
end


% Generate white noise
white = randn( pv.nrTimePoints,1);
white =(white-mean(white))/std(white);
% Fourier transform
[FT,freq] = fftReal(white,pv.samplingRate);
% Scale the amplitudes
scale = 1./(freq.^(pv.alpha/2)); %Divide alpha by two as alpha specifies drop off of power, and we're scaling amplitude
scale(1) =1; % Dont scale mean
FT = FT.*scale;
% Back to time domain
[signal,time] = ifftReal(FT,freq,pv.samplingRate);
% Force zero mean and unit stdev
signal = signal-mean(signal); % Zero center
signal = signal./std(signal); % sd =1

if nargout==0
    % Show
    subplot(1,2,1)
    plot(time,signal)
    xlabel 'Time'
    ylabel 'Signal'

    subplot(1,2,2)
    [FT,freq] = fftReal(signal,pv.samplingRate);
    
    plot(freq, FT.*conj(FT));
    set(gca,'YScale','log','XLim',[freq(2) freq(end)],'XScale','Log')
    title(['Power Spectrum of Colored Noise (\alpha = ', num2str(pv.alpha), ')']);
    xlabel('Frequency (Hz)');
    ylabel('Power');
end

end