function v = FBAR(v,tacsFrequency,samplingRate,frequencyLimits)
% FBAR - Remove tACS artifacts by removing the Fourier component 
% corresponding to the stimulation frequency (in each segment).
%
% INPUT
% v = segmented recorded signal [nrSamplesPerSegment nrSegments]
% tacsFrequency = Stimulation frequency in Hz.
% samplingRate = sampling rate of the recording in Hz.
% frequencyLimits = Frequencies within this band around the tacsFrequency
% will be removed. Defaults to [0 0] which means "the single frequency
% closest to the tacsFrequency". 
% 
% EXAMPLE tacsFrequency = 10 and frequencyLimits = [-1 1] would remove
% frequencies in the interval [9 11] Hz.
% OUTPUT
% v = segmented recovered signal.
arguments
    v (:,:) double 
    tacsFrequency (1,1) double
    samplingRate (1,1) double
    frequencyLimits (1,2) = [0 0];
end

% Vectorized FBAR
[ft, freq] = fftReal(v, samplingRate); % fftReal already handles matrices
if all(frequencyLimits==0)
    [~, ixTacs] = min(abs(freq - tacsFrequency));
else
    ixTacs = find(freq>=tacsFrequency+frequencyLimits(1) & freq<=tacsFrequency+frequencyLimits(2));
end
ft(ixTacs, :) = 0; % Zero out for all segments
v = ifftReal(ft, freq, samplingRate);

end
