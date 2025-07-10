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

for i=1:size(v,2)
    % For each segment ,do a FFT, find the relevant component and set it to
    % zero, then ifft.
    stay = ~isnan(v(:,i));
    [ft,freq] = fftReal(v(stay,i),samplingRate);
    [df,ixTacs] = min(abs(freq-tacsFrequency));
    if df ~=0
        fprintf(2,'No exactmatch for %.2f Hz tACS frequency. Removing %.2f Hz. Probably cause: non-integer number of cycles \n.',tacsFrequency,freq(ixTacs))
    end
    ft(ixTacs,:)=0;
    v(stay,i) = ifftReal(ft,freq,samplingRate);
end
