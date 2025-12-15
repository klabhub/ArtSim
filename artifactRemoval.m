function [vClean ,results] = artifactRemoval(v,pv)
% This function returns a cleaned up signal (vClean) following one or more
% artifact removal steps. 
%
% Briefly, there are the following steps (modes): (more details below)
% 'NOP' - no artifact removal. 
% 'NOTCH' - notch filtering
% 'RBAR' - regression based artifact removal
% 'FBAR'  - Fourier based artifact removal 
% 'FASTR' - Niazy et al algorithm as implemented in eegLab
% 'TC' - Time slice correction
% 'MEANREMOVAL' - remove the mean
%  'PCA'   - Use PCA to keep only some of the components.
% 'ANC-SEGMENT','ANC-ARTIFACT','ANC-REFERENCE' - Adaptive noise
%                   cancellation based on different models of the noise. 
% INPUT
% v = raw signal (assumed to be in volts) from which artifacts are to be removed.
% 'tacsFrequency' - Frequency of the tACS signal. It is assumed that tACS
% is applied for the duration of the signal v.
% 'segementDuration'   - For some artifact removal methods, the signal is divided
%               into segments.
% 'mode' - Cell array specifying which artifact removal algorithms to run
% (they will be executed in the specified order). For instance, the alogrithm
% in the Niazy et al paper (below) is approximated by: {'TC','MEANREMOVAL','PCA','ANC-ARTIFACT'}. 
% 
% Each mode has (partially overlapping) parameters, which are specified as parm/value
% pairs.
%% Mode: TC
%'Time Correction; If the event causing artifacts does not occur at
% regular times, Time Correction adjusts the segments used by MEANREMOVAL
% PCA, and ANC-ARTIFACT, such that they align best. This does not affect
% other ANC-modes.
%
% TC OPTIONS:
% 'slack' : Time in seconds over which the artifact onset time is expected
%                   to jitter [2e-3 s].
% 'referenceSegment': which segment to use as the reference. Defaults to 1
% (the first segment). If a vector of segments is specified, the average
% over those segments will be used as the reference. Specify 'mean' to use
% the mean over all segments as the reference.
%
%% Mode: 'MEANREMOVAL' 
% Remove the mean signal as determined in a window of segments
% around each segment.
% OPTIONS:
%  'nrSegPerWindow' - The window size in the number of segments. (1= one to
%                           the left and one to the right of the current segment)
% 'autocorrelationWindow' - Time in seconds over which autocorrelation in
%                           the signal is expected. Segments to determine the mean
%                           will be separated by at least this amount to avoid removing true signal.
%                           [0]
%  'slide'              Set to true to use a sliding window (width: nrSegPerWindow, step: 1 tACS cycle) to estimate the
%                           average artifact per segment. When false, the average artifact is determined as the average over all
%                           neighboring segments in the window. [true]
%% Mode : PCA
% OPTIONS:
%           pcaVarExplained - Keep the first N principal components that together explain
%                          at least this much variance. [75]
%           pcaNr         - Keep at most this number of components. [Inf].
%                           For FASTR mode only, set this to 'auto' to use the FASTR code
%                           to estimate the number of PCs to remove.
%           pcaZ           - Use only components with an explained variance
%                           that is an outlier (more than this number of standard deviations above the mean
%                           explained variance). [-inf]
%           pcaFastrUseSlice     - Only used in combination with FASTR mode; set to
%                           true to use slice-correction (instead of
%                           volume).
%           pcaNrSegsPerWindow' - Number of segments to use to determine an
%                   optimal basis set, using PCA. Larger windows give more accurate
%                   estimates, but at the cost of flexibility (i.e. they assume that the
%                   artifact is stationary on the scale of the window). (1= one to the left
%                   and one to the right). Defaults to Inf to include all segments to
%                   determine the basis set.
%
%  Note that the smallest number of PCs meeting each of the criteria is
%  kept.
%
%% Mode: ANC-ARTIFACT
%   Perform adaptive noise cancellation using the artifact signal
%   calculated up to that point (i.e. based on the modes preceding this
%   mode in the 'mode' argument).
% OPTIONS:
%   'ancMu' -  the step size of the ANC [0.05]
%   'ancOnline' - Whether to use online ANC filter estimates, or to
%   estimate the ANC filter based on the entire signal, then apply the
%   final estimate to all sapmples. [false].  (See anc.m for details)
%   For FASTR mode, setting ancOnline to false implies that ANC will not be
%   applied at all (it can only do the "online" variant).
%
%% Mode: ANC-SEGMENT
%  Perform adaptive noise cancellation using an artifact signal that
%  represents the onset of each segment (= a fixed phased of the tACS
%  signal).
% OPTIONS:
%   'ancMu','ancOnline'  - see above.
%% Mode: ANC-REFERENCE
%  With this mode the user can specif their own noise reference and use ANC to
%  remvove aspects of the signal that correlate with the noise reference.
% OPTIONS:
%  'ancReference', A vector representing the noise [nrSamples 1];
%   'ancN'  , the number of filter coefficients to use for ANC.
%   'ancMu','ancOnline'  - see above.
%   This mode can be added multiple times, with different noise references.
%   Note, however, that all ANC modes share ancMu and ancOnline, if
%   different mu/online parameters are needed, then call the artifactRemoval
%   function twice.
%% Mode: FASTR
%       This mode calls the fmrib eeglab plugin implementation of Nazy et
%       al's FASTR algorithm.  Note that there are a few idiosyncracies in
%       this implementation that may be ok for the intended application
%       (removing MRI artifacts) but not for the removal of tACS artifacts.
%       First, even when no PCA is requested (here: 'pcaExplainedVariance',0), 
%       the code regresses out the mean of each segment. This results in the 
%       removal of a lot of low frequency power (as is evident from the simulations in eegArtifacts.mlx).
%       Second, when PCA is requested (here:'pcaExplainedVariance','auto') it 
%       is performed on the high-pass filtered residuals (After removal of 
%        the mean artifact per segment).The filter cut-off is hard coded as 70 Hz.
%       Third, if the ANC is turned on (here: 'ancOnline',true),
%       the FASTR implementation automatically applies a 70 Hz
%       (hard-coded) low pass filter.
%% Mode: FBAR
% Fourier based artifact removal - removes all power that is exactly at the
% tACS frequency.
% OPTIONS:
% tacsFrequency - the frequency of the stimulation voltage
% 
%% Mode: RBAR
% Regression based artifact removal -  regresses out signal that is
% correlated with a specified voltage (the artifact proxy;  the measured
% voltage generated by the stimulator).
% OPTIONS
% vTacsRecord - the recorded tACS voltage.
% 
%% OTHER OPTIONS
% graph     -  show  a graph summarizing the artifact removal process.
%               [false]
% showTime   - Zoom in to the first few seconds in the graph [3]
% groundTruth - For simulated data, provide ground truth signal to
%               determine the quality of the artifact removal and report it on the
%               command line.
%
% OUTPUT
% vClean =  the signal after artifacts have been removed.
%
%
%% Sources:
% Allen, P. J., Josephs, O., & Turner, R. (2000). A method for removing imaging artifact from continuous EEG recorded during functional MRI.
% NeuroImage, 12(2), 230�239. https://doi.org/10.1006/nimg.2000.0599
%
% Niazy, R. K., Beckmann, C. F., Iannetti, G. D., Brady, J. M. & Smith, S. M.
% Removal of FMRI environment artifacts from EEG data using optimal basis sets. Neuroimage 28, 720�737 (2005).


arguments
    v (:,:) double
    pv.mode (1,:) cell {mustBeMember(pv.mode,{'NOP','NOTCH','RBAR','FBAR','FASTR','TC','MEANREMOVAL','PCA','ANC-SEGMENT','ANC-ARTIFACT','ANC-REFERENCE'})} = {}
    pv.segmentDuration (1,1) double = 3
    pv.tacsFrequency (1,1) double = NaN
    pv.autocorrelationWindow (1,1) double = 0
    pv.recordingSamplingRate (1,1) double = 30e3
    pv.slack (1,1) double = 2e-3
    pv.nrSegsPerWindow (1,1) double = 9
    pv.slide (1,1) logical = true
    pv.pcaVarExplained (1,1) double = 75
    pv.pcaNr = Inf
    pv.pcaZ (1,1) double = -inf
    pv.pcaFastrUseSlice (1,1) logical = false
    pv.pcaNrSegsPerWindow (1,1) double = inf
    pv.nrPCA (1,1) double = Inf
    pv.ancReference (:,1) double = []
    pv.ancN = []
    pv.ancOnline (1,1) logical = false
    pv.ancMu (1,1) double = 0.05
    pv.graph (1,1) logical = false
    pv.showTime (1,1) double = 2
    pv.groundTruth (:,:) double = []
    pv.vRecordTcs (:,:) double = []
    pv.referenceSegment = 1
end


doQC = ~isempty(pv.groundTruth);

nrSamples = size(v,1);
samplingRate= pv.recordingSamplingRate;
if pv.tacsFrequency>0 
    nrCycles = pv.segmentDuration*pv.tacsFrequency;
    if nrCycles~=round(nrCycles)
     fprintf(2,'Non-integer number of cycles per segment (%.2f).\n',nrCycles)
    end
else
    nrCycles = 1;
end
nrSamplesPerSegment =  pv.segmentDuration*samplingRate;
if nrSamplesPerSegment~=round(nrSamplesPerSegment)
    fprintf(2,'Non-integer number of samples per segment (%.2f).\n',nrSamplesPerSegment)
end
nrSegments = ceil(nrSamples/nrSamplesPerSegment);
if nrSegments~=nrSamples/nrSamplesPerSegment
    fprintf(2,'Non-integer number of segments (%.2f).\n',nrSamples/nrSamplesPerSegment);
end

% pad = zeros(samplesPerSegment,1);
% nrPad = numel(pad);
% v = [pad; v; pad];
nrPad=0;
segmentStart=nrPad+round(1+(0:nrSegments-1)*nrSamplesPerSegment);
[vSegmented,isPadded,nrPre,nrPost] = segment(v,segmentStart);


if ~isempty(pv.vRecordTcs)
    vTacsRecord = segment(pv.vRecordTcs,segmentStart);
else
    vTacsRecord = [];
end
if doQC
    groundTruth = pv.groundTruth;
end

vClean =vSegmented; % Start from the raw segmented.
totalArtifact = zeros(size(vClean));

if doQC
    [r,mse] = qc(v,groundTruth);
    if nargout >1
        results.raw.r = r;
        results.raw.mse = mse;
    else
        fprintf('********QC********* \n Starting Values: MSE = %3.3g muV , r= %3.3f\n',mse,r);
    end
end
modes = pv.mode;
nrModes = numel(modes);
filterToPlot = cell(1,nrModes);
artifactToPlot = cell(1,nrModes);

fourierFreq = pv.tacsFrequency;
for m=1:nrModes
    modeName = matlab.lang.makeValidName(modes{m});
    switch upper(modes{m})
        case 'NOP'
            % Do nothing
        case 'TC'
            % Do slice time correction to calculate artifact templates.
            if ~all(totalArtifact(:)==0)
                error('TC should be applied first.')
            end

            slackSamples = floor(pv.slack*samplingRate);
            if slackSamples > 0.5*ceil(samplingRate/pv.tacsFrequency)
                error('This much slack will lead to shifts larger than one cycle of tACS');
            end
            if strcmpi(pv.referenceSegment,'mean')
                reference = mean(vClean,2,'omitnan');
            else
                reference = mean(vClean(:,pv.referenceSegment),2,'omitnan');
            end
            unSegmentedVClean = unSegment(vClean,isPadded,nrPre,nrPost);
            shift = nan(1,nrSegments);
            for s =1:nrSegments
                % Pick samples around the (candidate) start of the segment
                pickIx = round((segmentStart(s)-slackSamples):(segmentStart(s)+nrSamplesPerSegment-slackSamples));
                pickIx(pickIx<1)=[];
                pickIx(pickIx>nrSamples)=[];
                candidate = unSegmentedVClean(pickIx);
                % Determine the best match with the reference
                [xc,lags] = xcorr(reference,candidate,slackSamples);
                [~,ix] = max(xc);
                nrPre = sum(pickIx<segmentStart(s));
                shift(s)= lags(ix)+nrPre;
            end
            segmentStart = segmentStart-shift; % Shift the segment start so that it matches best.
            [vClean,isPadded,nrPre,nrPost] = segment(unSegmentedVClean,segmentStart);
            totalArtifact = zeros(size(vClean));
            if nargout>1
                results.tcShift= shift;
            end
        case 'FBAR'
            vClean = FBAR(vClean,pv.tacsFrequency,samplingRate);
        case 'RBAR'
            vClean = RBAR(vClean,vTacsRecord);
        case 'NOTCH'
            d = fdesign.notch('N,F0,Q,Ast',500,pv.tacsFrequency,20,100,samplingRate);
            notch =design(d,'SystemObject',true);
            vClean = notch(unSegment(vClean,isPadded,nrPre,nrPost));
            vClean = segment(vClean,segmentStart);
        case 'MEANREMOVAL'
            % To include only segments where the underlying neural signal is assumed to
            % be uncorrelated; adjacent artifacts are not used but one or more are
            % skipped (step), depending on the assumed autocorrelationWindow
            step = 1+ceil(pv.autocorrelationWindow*samplingRate/nrSamplesPerSegment);
            halfWindow = pv.nrSegsPerWindow;
            meanArtifact = nan(size(vClean));
            removeMeanPerSeg = 0 ; % Set this to true to regress out the mean of each segment too. Not recommended (depending on the nrCycles used, this will also remove low frequency signals)
            scale = nan(1+removeMeanPerSeg,nrSegments);
            for s=1:nrSegments
                % Pick nrSegsPerWindow to estimate the artifact surrounding th current
                % segment (at the edges we only use the segments to one side of the segment that is to be corrected.)
                left = s - (step:step:2*halfWindow+1);
                right = s + (step:step:2*halfWindow+1);
                keep = [left(left>0) right(right<nrSegments)];
                [~,ix] =sort(abs(keep-s));
                keep  = keep(ix);
                tooMany = numel(keep)-(2*halfWindow+1);
                if tooMany>0
                    keep(end-tooMany+1:end)=[];
                end
                thisSegment    = vClean(:,s);
                stay = ~isnan(thisSegment);
                % Average over the artifact segments (vClean can have nans near the
                % edges)
                if pv.slide
                    % Sliding window
                    nrSamplesPerCycle = round(samplingRate/pv.tacsFrequency);
                    % Determine the mean artifact for each cycle within the
                    % window (i.e. ignoring the nrCycles that make up a segment)
                    singeCycleArtifact = mean(reshape(vClean(1:nrSamplesPerCycle*nrCycles,keep),nrSamplesPerCycle,[]),2,'omitnan');
                    % Then assume that this artifact simply repeats for all
                    % cycles within the segment.
                    ix = mod(0:numel(thisSegment)-1,nrSamplesPerCycle)+1;
                    meanTemplate   =  singeCycleArtifact(ix);
                else
                    meanTemplate = mean(vClean(:,keep),2,'omitnan');
                end
                % Determine how to scale this template to best match the current
                % segment
                dm = [ones(sum(stay),removeMeanPerSeg) meanTemplate(stay)];
                scale(:,s) =dm\thisSegment(stay);
                meanArtifact(stay,s) = dm*scale(:,s);
            end

            vClean = vClean- meanArtifact;
            totalArtifact = totalArtifact + meanArtifact;
            filterToPlot{m} = meanTemplate; %Show the last template
            artifactToPlot{m} = unSegment(meanArtifact,isPadded,nrPre,nrPost);

            if mean(abs(scale(end,:)-1))>0.05
                fprintf(2,'Large deviations from the template (%.2f%%)\n',100*mean(abs(scale-1)))
            end
            if nargout>1
                results.meanScale= scale;
                results.meanStdScale   = std(scale);
                results.meanMeanScale   = mean(scale);
            end

        case 'PCA'
            halfPCAWindow = min(pv.pcaNrSegsPerWindow,nrSegments);
            step = 1+ceil(pv.autocorrelationWindow/1000*samplingRate/nrSamplesPerSegment);
            pcaArtifact =zeros(size(vClean));
            cumVarExplained = zeros(1,nrSegments);
            nrPC = zeros(1,nrSegments);
            hasNan = find(any(isnan(vClean)));
            for s=1:nrSegments
                % Pick nrSegsPerWindow to estimate the artifacts surrounding the current
                % segment (at the edges we only use the segments to one side of the segment that is to be corrected.)
                left = (s-step):-step:(s-step*halfPCAWindow);
                right = (s+step):step:(s+step*halfPCAWindow);
                % For the PCA we randomly scatter some of the segments to include
                %They are always outside the autocorrelation window from the
                %current se
                if step>1
                    left = left-(randi(2,size(left))-1);
                    right = right + (randi(2,size(right))-1);
                end
                keep = [left(left>0) right(right<nrSegments)];
                keep = setdiff(keep,hasNan); % Remove segements with NaN (=padding)
                % Extract these artifacts surrounding the current segment
                artifactsForBasis =vClean(:,keep);
                if ~isempty(vTacsRecord)
                    artifactsForBasis = cat(2,artifactsForBasis,mean(vTacsRecord(:,keep),2,'omitnan'));
                end
                % Remove means from rows(samples) and columns (segments)
                % artifactsForBasis= artifactsForBasis-repmat(mean(artifactsForBasis),[nrSamplesPerSegment 1]) -repmat(mean(artifactsForBasis,2,'omitnan'),[1 numel(keep)]);
                % Perform PCA to find an optimal basis set to describe the current
                % segment
                [basis,~,~,~,varExplained] = pca(artifactsForBasis','Centered',false);
                % Selection
                lastPC = find(zscore(varExplained)>pv.pcaZ,1,'last');
                lastPC = min(lastPC,find(cumsum(varExplained)<pv.pcaVarExplained,1,'last'));
                lastPC = min(lastPC,pv.pcaNr);
                if ~isempty(lastPC)
                    cumVarExplained(s) = sum(varExplained(1:lastPC));
                    nrPC(s) = lastPC;
                    % Determine how to scale these PCs to best match the current
                    % segment (linear regression)
                    stay = ~(any(isnan(basis),2) | isnan(vClean(:,s)));
                    pcaScale = basis(stay,1:lastPC)\vClean(stay,s);
                    % Best esimate of the artifact
                    basis(~stay,1:lastPC) = 0;
                    pcaArtifact(:,s) = basis(:,1:lastPC)*pcaScale;
                end
            end
            % Remove the artifacts
            vClean = vClean- pcaArtifact;
            totalArtifact = totalArtifact + pcaArtifact;
            if ~isempty(basis)
                filterToPlot{m} = basis(:,1); % Show the first PC
            end
            artifactToPlot{m} = pcaArtifact(~isPadded);
            if nargout>1
                results.cumVarExplained = cumVarExplained;
                results.nrPC = nrPC;
            end
        case {'ANC-SEGMENT','ANC-ARTIFACT','ANC-REFERENCE'}
            % Do  ANC with a pulse at the start of each segment (i.e.
            % the start of a tACS sinusoid)
            switch upper(modes{m})
                case  'ANC-SEGMENT'
                    % Assume that each cycle of the tACS results in an
                    % identical artifact and remove that.
                    reference = zeros(nrSamples,1);
                    reference(segmentStart) =1;
                    % One coefficient for each time point in the segment
                    N = round(samplingRate.*(nrCycles/pv.tacsFrequency))-1;
                case 'ANC-ARTIFACT'
                    % Apply ANC to the artifact that has been estimated so
                    % far.
                    if all(totalArtifact==0)
                        N=0;
                    else
                        N = round(nrCycles.*samplingRate./pv.tacsFrequency)-1;
                        x = unSegment(vClean,isPadded,nrPre,nrPost);
                        y = unSegment(totalArtifact,isPadded,nrPre,nrPost);
                        reference = ((x'*y)./(y'*y))*y;
                    end
                case 'ANC-REFERENCE'
                    % This ANC variant uses a reference [nrSamples 1] provided by the caller. In principle, the caller
                    %  can specify any reference here to remove dataq that correlate with this reference from the signal.
                    N =   pv.ancN;
                    reference =pv.ancReference./max(abs(pv.ancReference));
            end
            if N>0
                [vClean,ancNoise,ancFilter]=anc(reference,unSegment(vClean,isPadded,nrPre,nrPost),N,pv.ancMu,[],pv.ancOnline);
                totalArtifact = totalArtifact+ segment(ancNoise,segmentStart);
                vClean= segment(vClean,segmentStart);
                filterToPlot{m} = ancFilter;
                artifactToPlot{m} = ancNoise;
            else
                filterToPlot{m} = 0;
                artifactToPlot{m} = zeros(nrSamples,1);
            end
        case 'FASTR'
            % FMRI artifact Slice Template Removal (FASTR) as implemented
            % in fmrib eeglab plugin.
            % From fmrib_fastr:
            %   EEG:  EEGLAB data structure
            %   lpf:  low-pass filter cutoff
            %   L: Interpolation folds
            %   window: length of averaging window in number of artifacts
            %   Trigs: An array of slice triggers locations.
            %   strig: 1 for slice triggers, 0 for volume / section triggers.
            %   anc_chk: 1 to do Adaptive noise cancellation
            %			 0 to not.
            %   tc_chk:  1 to correct for missing triggers, 0 for not
            %   Volumes: FMRI volumes for use in trigger correction
            %   Slices:  FMRI Slices / Volume for use in trigger correction
            %   varargin{1}: relative position of slice trigger from beginning of
            %       slice acquisition: 0 for exact start -> 1 for exact end
            %       default=0.03;
            %   varargin{2}: Channels not to perform OBS  on.
            %   varargin{3}: Numer of PCs to use in OBS. use 0 to skip this step.
            %                'auto' or empty for auto order selection.
            %

            %% Re-package parameters to match fmriB code
            EEG.data = v';           % The data [1 nrSamples]  (1 = nr Channels)
            EEG.srate = samplingRate;  % Sampling rate
            EEG.pnts= size(EEG.data,2);            % Number of data points
            % Low-pass filter cut-off. Applied after mean removal and PCA, before ANC.
            % With this set to 0 (no filter) ANC will take a long time, but
            % at least we see the full result.
            lpf = 0;
            L = 1;  % "Interpolation folds" = upsampling. We are already at 10kHz so not needed.
            Window = 2*pv.nrSegsPerWindow; % Number of segments to average for one template.            
            % The realignment algorithm needs time before the first and after the last
            % trigger to work.
            pad = zeros(1,nrSamplesPerSegment); %
            nrPad = numel(pad);
            Trigs =nrPad+(1:nrSamplesPerSegment:nrSamples);
            EEG.data = [pad EEG.data pad];
            strig = pv.pcaFastrUseSlice; % 1= slice triggers (0 =volume triggers)
            anc_chk = pv.ancOnline; % Do ANC or not
            tc_chk = 0;  % Option to correct missing triggers. Not relevant
            Volumes =[];  % Only used to correct missing triggers
            Slices = [];  % Only used to correct missing triggers
            % The realignment slack is coded as a fraction of the time
            % between triggers:
            pre_frac =pv.slack/(nrSamplesPerSegment/samplingRate);
            noObsChannels =[];  % Channels to exclude from PCA (none here)
            % 0 = no PCA, 'auto' = determine number of PC based on data, n= use this number of components.
            nrPCA =pv.pcaNr;

            %% Now call the original code
            results = fmrib_fastr(EEG,lpf,L,Window,Trigs,strig,anc_chk,tc_chk,Volumes,Slices,pre_frac,noObsChannels,nrPCA);
            [vClean,isPadded,nrPre,nrPost] = segment(results.data(nrPad+1:end-nrPad)',segmentStart);
            filterToPlot{m} = 0;
            artifactToPlot{m} = zeros(nrSamples,1);
            %case 'MYMODE'
            % Extend ArtSim with novel methods by adding artifact removal
            % steps here . Note that the arguments block needs to be told about
            % the mode too (see top of function).
        otherwise
            error('Unknown artifact removal mode %s',modes{m})
    end

    if doQC
        % Show quality of the reconstruction (mse and correlation)
        [r,mse] = qc(unSegment(vClean,isPadded,nrPre,nrPost),groundTruth,fourierFreq,samplingRate);
        if nargout >1
            results.r.(modeName) = r;
            results.mse.(modeName) = mse;
        else
            fprintf('\t After %s: MSE = %3.3g muV , r= %3.3f\n',modes{m},mse,r);
        end
    end

end

if nargout>1
    results.modes= modes;
    results.artifacts = artifactToPlot;
    results.filters= filterToPlot;
end

%Reshape to the original size.
vClean = unSegment(vClean,isPadded,nrPre,nrPost);
assert(numel(vClean)==nrSamples,'Losing samples...');

%% Show results 
if pv.graph
    plotModes= find(~cellfun(@isempty,artifactToPlot));
    clf;
    t = (0:nrSamples-1)/samplingRate;

    %% Artifacts detected by each mode
    subplot(1,3,1);
    hold on
    for m=plotModes
        plot(t,artifactToPlot{m}/1e-3);
    end
    xlabel 'Time (s)'
    ylabel 'Voltage (mV)'
    xlim([0 pv.showTime]);
    title 'Artifact (mV)'
    legend(modes(plotModes))

    %% Filters created in each mode
    subplot(1,3,2)
    hold on
    norm = @(x) (x./max(abs(x)));
    nrm = nan(1,numel(plotModes));
    for m=plotModes
        plot(norm(filterToPlot{m}))
        nrm(m) = max(abs(filterToPlot{m}));
    end
    xlabel 'Samples'
    ylabel 'Relative magnitude'
    title(['Artifact Filters/Templates' num2str(nrm,2)])
    legend(modes(plotModes))

    % Raw and cleaned signal
    subplot(1,3,3)
    plot(t,v(:)/1e-3,'LineWidth',1)
    hold on
    plot(t,vClean/1e-3,'LineWidth',2);
    if doQC
        plot(t,groundTruth(:)/1e-3,'g','Linewidth',1);
    end
    xlabel 'Time (s)'
    ylabel 'Voltage (mV)'
    xlim([0 pv.showTime]);
    title 'Recorded and cleaned signal'
end
end

%% Helper functions
function v=unSegment(v,isPadded,nrPre,nrPost)
v=[zeros(nrPre,1);v(~isPadded);zeros(nrPost,1)];
end

function [v,isPadded,nrPre,nrPost] =segment(v,segmentStart)
nrSamples= numel(v);
[sampleIx,isPadded,nrPre,nrPost] = splitInSegments(nrSamples,segmentStart);
v =[v;NaN];
v=v(sampleIx);
end

function [sampleIx,isPadded,nrPre,nrPost] = splitInSegments(nrSamples,segmentStart)
if isscalar(segmentStart) && segmentStart==1
    maxNrSamplesPerSegment = nrSamples;
    nrSamplesPerSegment = nrSamples;
else
    maxNrSamplesPerSegment = max(diff(segmentStart));
    nrSamplesPerSegment = [diff(segmentStart) maxNrSamplesPerSegment];
end
tmp = arrayfun(@(x) (x+(0:maxNrSamplesPerSegment-1)'),segmentStart,'Uni',false);
isSpillOver = arrayfun(@(x) ([false(x,1);true(maxNrSamplesPerSegment-x,1)]),nrSamplesPerSegment,'Uni',false);
isSpillOver = cat(2,isSpillOver{:});
sampleIx= cat(2,tmp{:});

lastIxInUse = max(sampleIx(:));
sampleIx(sampleIx>nrSamples) = nrSamples+1;% Beyond last sample
sampleIx(sampleIx<1) = nrSamples+1; % Before first sample
isPadded= sampleIx ==nrSamples+1 | isSpillOver;
nrPre = max(0,segmentStart(1)-1);
nrPost = max(0,nrSamples-lastIxInUse);
end