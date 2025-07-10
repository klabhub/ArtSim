function [filtered,noiseEstimate,filter]= anc(reference,signal,N,mu,guess,online)
% Adaptive Noise Cancellation
%
% INPUT
% reference -  the reference signal (i.e. some signal that represents noise
%               to be removed).  [nrSamples 1]
% signal    - the signal from which noise is to be removed [nrSamples 1]
% N         - the number of filter coefficients to use [20]
%
% mu        - The "step size", that represents by how much the filter is
%               updated each step. Note that this number will be scaled by
%               1./(N*stdev(reference)) using the MAD estimate for the
%               stdev. [0.05]
% guess      - The initial guess for the filter (set to [] to start with
%              zeros).
% online      - Set to false to apply the final filter estimate to all data samples.
%               [true]
% OUTPUT
% filtered   - The signal, after removing the noise as estimated by ANC
% noise      - The estimated noise
% filter     - The estimated filter
% delta       - the mean absolute change in the filter with each time step.
%
% BK  - Nov 2021

nin=nargin;
if nin <6
    online = true;
    if nin <5
        guess = [];
        if nin<4
            mu =0.05;
            if nin <3
                N =20;
            end
        end
    end
end

%#ok<*ASGLU> 
%#ok<*UNRCH> 

debug =false; 
useMex = false; % The Matlab code below
                % is faster than the compiled mex file...so don't
                % use mex...

assert(size(signal,2)==1 && size(reference,2)==1,'Only single columnvectors are allowed in the call to anc.m')
assert(size(signal,2)==size(reference,2),'Reference and Signal must have the same number of samples')
assert(~any(isnan(reference)),'The reference signal should not have NaNs.')
nrSamples=numel(reference);
% Scale mu as in Allen et al 2010 but with a robust estimate of the stdev
mu = mu./(N.*mad(reference));
if isempty(guess)
    guess = zeros(N+1,1);
end

if  useMex && exist('ancEstimate','file')==3
    % Use the mex-file to estimate the filter. Same results.
    [onlineFiltered,noiseEstimate,filter,delta] =ancEstimate(reference,signal,N,mu,guess);    
else    
    % Use the same algorithm in Matlab.
    if debug
        everyN =100; 
        allFilters= nan(N+1,ceil(nrSamples/everyN));
        cntr= 0;
    end
    %Initialize
    thisReference=flipud([0;reference(1:N)]);
    filter=guess;
    noiseEstimate =zeros(nrSamples,1);
    onlineFiltered = signal;
    delta = nan(nrSamples,1);
    % Estimate the optimal filter by sliding over the signal & reference
    for sample =N+1:nrSamples
        thisReference=[reference(sample); thisReference(1:end-1)];
        noiseEstimate(sample)=filter'*thisReference;
        onlineFiltered(sample) =signal(sample)-noiseEstimate(sample);
        % Update rule
        filterChange = 2*(onlineFiltered(sample)*mu).*thisReference;        
        delta(sample) = mean((filterChange));
        filter = filter + filterChange;
        if debug 
            if mod(sample,everyN)==0
                cntr= cntr+1;
                allFilters(:,cntr) = filter;
            end
        end
    end    
end

if online
    % Use the online filtering (as estimated above)
    % This makes sense if the artifact is non-stationary on the timescale of
    % the signal provided here. The cost is that the start of the signal
    % segment will inevitably have a poor esitimate/noise correction
    % (unless convergence is really rapid, or the initial guess is close to
    % the true filter).
    filtered = onlineFiltered;
else
    % Apply the final filter estimate to all of the refernce and subtract from
    % the original signal.
    % This makes most sense if the artifact is stationary and the only
    % variation in the filter is due to the estimation process. The advantage 
    % is that even the earliest samples will be corrected. 
    singleFilterNoiseEstimate = conv(reference,filter,"full");
    singleFilterNoiseEstimate(end-N+1:end) = [];
    filtered  =signal - singleFilterNoiseEstimate;
end
end