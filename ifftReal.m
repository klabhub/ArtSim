function [signal,time] = ifftReal(ft,freq,sf,N)
% Inverse Fourier of the lower half of the spectrum (i.e. the inverse of
% fftReal). The ft input can be a matrix, but the first dimension has to
% contain the components for the different frequencies.
% 
% INPUT
% ft - Fourier components as returned by fftReal
% freq - Frequencies of the ft (as returend by fftReal)
% sf  - sampling frequency
% N = the number of samples to use in the actual call to IFFT (see help IFFT
% for more information about the truncation or zero-padding that is used). Defaults to
% [] (i.e. use all samples, wihtout zero-padding)
%
% OUTPUT
% signal - The signal in the time domain
% time - Corresponding time points.

if nargin<4
    N = [];
    if nargin <3 
        sf= 1;
    end
end

nrComps =size(ft,1);
if numel(freq) ~=nrComps
    error('Frequencies should be along row dimension');
end
% Let's make sure the freq and ft are sorted correctly.
[freq,ix] =sort(freq,'ascend');
ft =ft(ix,:);
   
% Determine whether the highest frequency is the Nyquist (even samples in
% the original signal) or a frequency below (odd samples).
hasNyquist = abs(freq(end)-sf/2)<eps; % If the difference is smaller than eps, they are really the same to machine precision.

% Determine the "mirrored" components using the complex conjugate
if hasNyquist
    % The last component is the
    % Nyquist, don't mirror that one. 
    % Don't mirror the DC
    mirror = conj(ft(end-1:-1:2,:));      
    nrSamples = nrComps+size(mirror,1); 
    % From fftReal we have :
    % highestFreq = nrSamples/2*sf/nrSamples;  
    % highestFreq =  sf/2   
    % Therefore 
    computedSamplingFrequency = 2*freq(end);      
    % (which we could also have derived from sf = 2*nyquist ).
else
    % The last comp is below nyquist, mirror that one.
    % Don't mirror the DC.
    mirror = conj(ft(end:-1:2,:));
    % From fftReal we have
    % highestFreq = ((nrSamples-1)/2)*sf/nrSamples;
    % Inverting this :
    nrSamples = nrComps+size(mirror,1);
    computedSamplingFrequency = freq(end)*nrSamples/((nrSamples-1)/2);    
end
% Combine the provided Fourier
% components, and the mirrored components to make the fullFt symmetric
fullFt = cat(1,ft,mirror);
% Now we can use ifft along dimension 1
signal = ifft(fullFt,N,1);

deltaSf = abs(sf-computedSamplingFrequency);
if deltaSf>eps
    warning('The specified sampling frequency (%3.2g) and the computed sampling frequency (%3.2g) are mismatched',sf,computedSamplingFrequency);
end
dt = 1./sf;
time = (0:nrSamples-1)*dt;

end

