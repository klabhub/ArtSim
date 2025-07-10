function [ppcValue,pValue] = ppc(spikeIx,phase,nrBootstrap)
% function [ppcValue,pValue] = ppc(spikeIx,phase,nrBootstrap)
% Estimate the pairwise phase consistency as defined by Vinck et al. 
% The Phase Locking Value (PLV) is the square root of the PPC.
%
% Note that a parfor is used for bootstrapping the p-value; for nrBoostrap>0
% this function will start the default parallel workers on your machine.
%
% INPUT
% spikeIx - Index into the phase vector where spikes occurred. This can be
% a cell array with two such vectors, then the PPC is computed for each,
% and the bootstrap test assesses the significance of their difference.
% phase   - phase at each sample point
% nrBootstrap - Number of boostrap samples to use for significance testing.
% [0]
%  OUTPUT
% ppcValue - The PPC (The PLV is the sqrt of this).
% pValue - The corresponding p-value

arguments
    spikeIx
    phase (:,1) double
    nrBootstrap (1,1) double = 0
end
if iscell(spikeIx) 
    % Evaluate the PPC/PLV for the two vectors of spike indices in the cell array,
    % then determine whether these PLVs are significantly different.     
    ppcStim = locPpc(spikeIx{1},phase(:,1));
    ppcNoStim = locPpc(spikeIx{2},phase(:,2));
    deltaPlv = abs(ppcStim-ppcNoStim);
    ppcValue = [ppcStim ppcNoStim];

    if nargout>1  && nrBootstrap>0
        % Bootstrap for difference
        allSpikeIx = [spikeIx{1} spikeIx{2}];
        nrNoStim = numel(spikeIx{2});
        nrStim = numel(spikeIx{1});
        nrAll = nrStim+nrNoStim;
        deltaBs = nan(1,nrBootstrap);
        parfor i=1:nrBootstrap
            a = locPpc(allSpikeIx(randperm(nrAll,nrStim)),phase);
            b = locPpc(allSpikeIx(randperm(nrAll,nrNoStim)),phase);
            deltaBs(i) = abs(a-b);
        end
        pValue = max(1/nrBootstrap,1-sum(deltaPlv>deltaBs)/nrBootstrap);
    else
        pValue = NaN;
    end
else
    % Asses PPC for a single set of spike times and determine whether this
    % value is significantly different from zero.
    ppcValue = locPpc(spikeIx,phase);
    if nargout>1  && nrBootstrap>0
        % Bootstrap for difference from random
        bsPpc = nan(1,nrBootstrap);
        nrSpikes= numel(spikeIx);
        parfor i=1:nrBootstrap
            randIx = randperm(numel(phase),nrSpikes);
            bsPpc(i) = locPpc(randIx,phase);                      
        end
        pValue = max(1/nrBootstrap,1-sum(ppcValue>bsPpc)/nrBootstrap);
    else
        pValue = NaN;
    end
end

end

%% Function that computes PPC.
function v = locPpc(spikeIx,phase)
nrSpk = numel(spikeIx);
if nrSpk>1
    spikePhase =  phase(spikeIx);
    cPhase = cos(spikePhase)';
    sPhase = sin(spikePhase)';
    tmp = 0;
    for i=1:nrSpk-1
        for j= ((i+1):nrSpk)
            tmp = tmp + cPhase(i)*cPhase(j) + sPhase(i)*sPhase(j);
        end
    end
    v = (2/(nrSpk*(nrSpk-1))*tmp);
else
    v =0;
end
end

