function v = FBAR(v,tacsFrequency,samplingRate)
% FBAR - Remove tACS artifacts by removing the Fourier component 
% corresponding to the stimulation frequency (in each segment).
%
% INPUT
% v = segmented recorded signal [nrSamplesPerSegment nrSegments]
% tacsFrequency = Stimulation frequency in Hz.
% samplingRate = sampling rate of the recording in Hz.
% OUTPUT
% v = segmented recovered signal.
arguments
    v (:,:) double 
    tacsFrequency (1,1) double
    samplingRate (1,1) double
end

% Vectorized FBAR
[ft, freq] = fftReal(v, samplingRate); % fftReal already handles matrices
[~, ixTacs] = min(abs(freq - tacsFrequency));
ft(ixTacs, :) = 0; % Zero out for all segments
v = ifftReal(ft, freq, samplingRate);
end
