function [ft,freq] = fftReal(signal,sf,N)
% function [ft,freq] = fftReal(signal,sf,N); 
% Return only the first half of the Fourier components (and their
% frequencies). The samples are assumed to be in dimension 1 of signal. 
%
% INPUT
% signal =  the signal
% sf = sampling frequency. Defaults to 1.
% N = the number of samples to use in the actual call to FFT (see help FFT
% for more information about the truncation or zero-padding that is used). Defaults to
% [] (i.e. use all samples, wihtout zero-padding)
% OUTPUT
% ft  = Fourier components
% freq = Corresponding frequencies
%
% NOTE
% Use ifftReal to transform these ft back to the time domain.

if nargin<3
    N=[];
    if nargin <2
        sf = 1;
    end
end
if any(~isreal(signal))
    error('fftReal only works for signals without complex values');
end

nrSamples =size(signal,1);
if rem(nrSamples,2)==0
    % Even number of samples. Nyquist is the highest one
    frequency = ([0 1:nrSamples/2 (nrSamples/2-1):-1:1]*sf/nrSamples)';  
    [~,highestFrequencyIx] =max(frequency);
else
    % Odd number of samples. Nyquist is not in the set.
    frequency = ([0 1:(nrSamples-1)/2 (nrSamples-1)/2:-1:1]*sf/nrSamples)';
    [~,highestFrequencyIx] =max(frequency);
end
% Select the lower half
keep = 1:highestFrequencyIx;
freq = frequency(keep);
signal(isnan(signal))=0;

% Force fft to work along dimension one of the signal matrix.
ft = fft(signal,N,1);

% This is a trick that allows fftReal to work along columns of a matrix
% or higer dimensional arrays, essentialy we're constructing a line like
% ft = ft(keep,:,:). 
allDims(1:ndims(signal)-1) = {':'};
ft = ft(keep,allDims{:});