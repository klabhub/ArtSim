classdef artSim < handle
    % artSim  A Matlab handle class to simulate electrophysiological
    % recordings, investigate the influence of artifacts induced by electrical
    % stimulation, and test artifact removal algorithms.
    %
    % By default class members are set to have "no influence". E.g, a
    % a spike/lfp/tacs amplitude of 0, a 1000V linear amplification
    % range,etc.
    properties (SetAccess=public,GetAccess=public)
        reproducible              = false;    % Set to true to reset the RNG to default before each simulation.
        % Set it to a positive integer to use that as a seed.

        duration                  = 120;       % [s] Duration of the recording to simulate.
        fromData                  = false;     % If lfp/eeg data are extracted from a data file, this is set to true to skip simulation.

        %% Stimulation
        recordingSamplingRate     = 30e3;     % [Hz] Sampling rate of the recording device
        stimulatorSamplingRate    = 2e3;      % [Hz] Sampling rate of the TCS stimulator
        currentResolution         = 1e-6;     % [A] Current resolution of the stimulator.

        % tACS
        tacsAmplitude             = 0;        % [A] Current amplitude
        tacsPhaseOffset           = 0;        % [Rad] Current phase
        tacsFrequency             = 0;        % [Hz] Current frequency
        tacsShape                 = 'sine';   % ['sine'] {'sine','sawtooth','square','none','am'}
        tacsCarrierFrequency      = 0;        % [Hz] Carrier frequency for AM Modulated tACS.
        % tDCS
        tdcsMean                  = 0;       % [A] Current
        % Ramp (applied to both tACS and tDCS)
        tcsRamp                  = 0;       % [s] Duration of the on and offset ramp.

        resistance                = 1;        % [Ohm] Effective resistance at the recording site.

        %% Physiological impedance artifacts
        % Simulating periodic resistance changes as caused by physiological changes such as heartbeats or respiration
        % The default implementation uses a periodic sinusoidal or
        % pulse-like shape, but the user can specify any zFun to simulate
        % more complex shapes (see zHeartbeat.m for an example)
        zFrequency        = 0;        % [Hz]
        zAmplitude        = 0;        % [fraction of the resistance]
        zDuration         = 0.25      % [s] - Full width at half max.
        zVariability      = 0;        % [fraction of heartbeat frequency per *hour*]
        zFun             = []         % Function that returns the impedance at each timepoint.

        %% Spikes
        spikePeak                 = 0;        % [V] Peak of the action potential
        spikePeakVariation        = 0;        % Spikes are scaled by 1+randn*spikePeakVariation.
        spikeRate                 = 0;        % [Hz] Firing rate
        poisson                   = true;     % Set to true to generate spike times with a Poisson distribution.
        % When false, spike times are drawn from a flat distribution.
        fixNrSpikes               = false;    % Fix the number of spikes to match the rate. (see simSpikes)
        refractoryPeriod          = 1.5e-3;   % [s] Refractory period
        synchrony                 = 0;        % Von Mises Kappa:
        % 0 means Poisson random spike times,
        % Typical PLV <0.05 which corresponds to a synchrony of 0.1
        % A high PLV of 0.3 is reached for synchrony=1;
        %% Local Field Potentials and EEG signals
        lfpAmplitude              = 0;        % [V] amplitude of the intrinsic LFP signal.
        lfpPhaseOffset            = 0;        % [rad]
        lfpFrequency              = 0;        % [Hz] Frequency of an intrinsic oscillation
        lfpPeakWidth              = 0;        % [Hz] Standard deviation of the amplitude distribution for a collection of LFPs with random phases.
        preferredPhase            = 0;        % [Rad] Phase at which spikes occur preferentially (if synchrony>0)

        %% Background activity/noise
        additive                  = 0;        % [V]  Standard deviation of additive noise
        pinkNoise                 = true;     % Set to true to generate pink noise (white otherwise)
        multiplicative            = 0;        % Standard deviation of the multiplicative noise

        %% Properties of the amplifier and ADC converter.
        % amplifierLinearRange is the output-referred linear range of the
        % amplifier. The input-referred range is amplifierLinearRange /
        % amplifierGain.
        amplifierGain             = 1;         % Preamp gain.
        amplifierLinearRange      = 100e-3;    % [V] Output-referred linear range
        adcNonlinearity           = 100;       % The saturation profile of the amplifier near the edge
        %                                       of the linear range. 2 = very sigmoidal, 100 = linear inside the range,
        %                                       saturation outside (i.e. piecewise linear). See amplify
        adcRange                  = 100e-3;     % The range of the ADC; bits are distributed evenly across this +/- range.
        adcBitDepth               = 24;       % Number of bits in the ADC


        highPass                  = 0;        % [Hz] Analog highpass cutoff before amplification.
        lowPass                   = 0;        % [Hz] Analog anti-alias lowpass cutoff before the ADC.
                                              % If this is zero, a built-in antialiasing filter is used
                                              % before downsampling to the recording sampling rate.
        digitalLowPass            = 0;        % [Hz] Optional lowpass cutoff after ADC quantization.
    end
    properties (SetAccess= protected, GetAccess=public)
        tSimulate;  % The internal time of simulation
        ixSpike;    % Index (into vNeural/tSimulate) where spikes occur
        vSpike;     % The extracellular voltage produced by the spikes
        vSpikeNeighbor; % The extracellular voltage produced by (other) spikes at a neighboring electrode.
        lfp;        % The LFP at the electrode of interest.
        vTcsActual % The voltage created by the TCS (includes discretization and impedance artifacts)
        tRecord;     % The time of each sample
        vRecord;    % The recorded voltage - including the artifact. First column is without artifact , second column with.
        vRecordNeighbor; % The recorded voltage at the "neighboring" electrode.First column is without artifact, second column with.
        vRecordTcs; % The recorded TCS voltage (at the recording sampling rate).
    end

    properties (Dependent)
        tTcs       % The time of each stimulator sample
        nrSamples;  % The number of simulated recorded samples (at recordingSamplingRate)
        nrTimePoints; % The number of simulated time points ( at quasiContinuousSamplingRate)
        phaseLfp;   % The phase of the LFP at each time point.
        vNeural;    % The neural voltage trace (spikes + lfp)
        vNeuralNeighbor   % The neural voltage trace at a neighboring electrode.
        vTruth;     % The true voltage (as would have been recorded without the artifact)
        vContaminated; % The recorded voltage (with the artifact).
        quasiContinuousSamplingRate;  % "Sampling" rate of the underlying simulation
        phaseTacs; % The phase of the recorded tACS signal
        nrSpikes;  % Total number of simulated spikes
        adcVoltsPerBit; % Volts for each bit in the ADC
    end

    methods  %get/set for dependent properties
        function set.duration(o,v)
            o.duration = v;
            reset(o)
        end
        function set.recordingSamplingRate(o,v)
            o.recordingSamplingRate = v;
            reset(o)
        end
        function v = get.quasiContinuousSamplingRate(o)
            % Integer multiple of the recording sampling rate, above 100
            % kHz to avoid sharp simulation edges to show up in recorded
            % data.
            overSampling = ceil(100e3/o.recordingSamplingRate);
            v = overSampling*o.recordingSamplingRate;
        end

        function v = get.phaseTacs(o)
            v = artSim.phase(o.vRecordTcs);
        end
        function v = get.nrSpikes(o)
            v= numel(o.ixSpike);
        end
        function v= get.vTruth(o)
            % Returns the recorded voltage as would have been recorded
            % in the absence of artifacts.
            v = o.vRecord(:,1);
        end

        function v= get.vContaminated(o)
            % Returns the recorded voltage including all artifacts
            v = o.vRecord(:,2);
        end
        function v= get.tTcs(o)
            % Returns the simulated time points (at the stimulator sampling
            % rate)
            v = (0:1/o.stimulatorSamplingRate:(o.duration-1/o.stimulatorSamplingRate))';
        end
        function v= get.nrSamples(o)
            % Returns the number of samples
            v = o.duration*o.recordingSamplingRate;
        end
        function v= get.nrTimePoints(o)
            % Returns the number of time points in the simulation
            v = o.duration*o.quasiContinuousSamplingRate;
        end
        function v= get.vNeural(o)
            % Returns the voltage at the electrode of interest.
            v = o.vSpike + o.lfp;
        end
        function v= get.vNeuralNeighbor(o)
            % Returns the  voltage at a neighboring electrode
            v = o.vSpikeNeighbor+o.lfp;
        end

        function v= get.phaseLfp(o)
            % Returns the phase of the LFP for each time point.
            v= artSim.phase(o.lfp);
        end

        function v = get.adcVoltsPerBit(o)
            % The size of each bit step in volts.
            v = 2*o.adcRange/(2^o.adcBitDepth-1);
        end
    end

    methods (Access=public)
        function disp(o)
            %disp - Display public simulation parameters by concept
            fprintf('artSim simulation parameters\n');
            displayGroup('Simulation', ...
                {'reproducible','duration','fromData'});
            displayGroup('Recording and stimulation', ...
                {'recordingSamplingRate','stimulatorSamplingRate', ...
                'currentResolution','resistance','tacsAmplitude', ...
                'tacsPhaseOffset','tacsFrequency','tacsShape', ...
                'tacsCarrierFrequency','tdcsMean','tcsRamp'});
            displayGroup('Physiological impedance artifacts', ...
                {'zFrequency','zAmplitude','zDuration','zVariability','zFun'});
            displayGroup('Spikes', ...
                {'spikePeak','spikePeakVariation','spikeRate','poisson', ...
                'fixNrSpikes','refractoryPeriod','synchrony','preferredPhase'});
            displayGroup('Local field potential', ...
                {'lfpAmplitude','lfpPhaseOffset','lfpFrequency','lfpPeakWidth'});
            displayGroup('Background noise', ...
                {'additive','pinkNoise','multiplicative'});
            displayGroup('ADC and recording filters', ...
                {'amplifierLinearRange','adcNonlinearity','adcRange', ...
                'adcBitDepth','amplifierGain','highPass','lowPass', ...
                'digitalLowPass'});

            function displayGroup(title,names)
                fprintf('\n%s:\n',title);
                for k = 1:numel(names)
                    displayProperty(names{k});
                end
            end

            function displayProperty(name)
                value = o.(name);
                if isa(value,'function_handle')
                    valueText = ['@' func2str(value)];
                elseif ischar(value)
                    valueText = ['''' value ''''];
                elseif isstring(value)
                    valueText = ['"' char(value) '"'];
                elseif isnumeric(value) || islogical(value)
                    valueText = mat2str(value);
                else
                    valueText = class(value);
                end
                fprintf('  %-24s %s\n',[name ':'],valueText);
            end
        end

        function reset(o)
            % Reset the object.
            o.ixSpike = [];
            o.vSpike = zeros(o.nrTimePoints,1);
            o.vSpikeNeighbor = zeros(o.nrTimePoints,1);
            if ~o.fromData
                o.lfp = zeros(o.nrTimePoints,1);
            end
            o.tRecord = (0:1/o.recordingSamplingRate:(o.duration-1/o.recordingSamplingRate))';
            o.tSimulate = (0:1/o.quasiContinuousSamplingRate:(o.duration-1/o.quasiContinuousSamplingRate))';
            o.vRecord = zeros(o.nrSamples,2);
            o.vRecordNeighbor = zeros(o.nrSamples,2);
            o.vRecordTcs = zeros(o.nrSamples,1);
        end
        function o = artSim(prms)
            %artSim Construct an artSim instance
            % Can be called with a structure that has fields matching
            % artSim properties to set those properties.
            arguments
                prms (1,1) struct = struct; % Parameter settings.
            end
            if ~isempty(prms)
                % Must be a struct with field names that match properties
                % of the artSim class.
                fn = fieldnames(prms);
                for p=1:numel(fn)
                    o.(fn{p}) = prms.(fn{p});
                end
            end
            reset(o);
        end

        function rng(o)
            % This is called by the simPlot functions to reset the RNG.
            % With .reproducible= true, the same results are regenerated
            % (useful for debugging or reporting results in a paper).
            if islogical(o.reproducible)
                if o.reproducible
                    rng default
                else
                    % false - do nothing
                end
            elseif isnumeric(o.reproducible) && o.reproducible >0 && round(o.reproducible)==o.reproducible
                % Use as seed
                rng(o.reproducible)
            else
                error('reproducible should be true (reset rng to default), false (keep rn as is) or a positive integer (the seed)')
            end
        end

        function data = addData(o,src,pv)
            % This function loads EEG or LFP data to serve as ground truth.
            % The src  can be an EEGLAB EEG struct or a Fieldtrip data struct
            % from which the first o.duration seconds will be used.
            % Alternatively, you can supply a vector of samples that
            % matches o.duration*o.recordingSamplingRate ("RAW").
            %
            arguments
                o (1,1) artSim
                src  % The EEGLAB EEG struct or a Fieldtrip data struct, or a matrix of doubles with samples along the first dimension, channels along the second
                pv.format (1,1) string {mustBeMember(pv.format,["EEGLAB" "FIELDTRIP" "RAW"])} = "EEGLAB";
                pv.channel (1,:) =1 % Which channel(s) to use from a multichannel dataset (will be averaged)
            end
            switch pv.format
                case "EEGLAB"
                    % Read the first o.duration  seconds from continuous
                    % data stored in an EEGLAB struct
                    o.recordingSamplingRate = src.srate;
                    assert(o.duration<=src.xmax,"Not enough samples in this EEGLAB struct");
                    src = pop_select(src,'time',[0 o.duration-1/src.srate],'channel',pv.channel);
                    o.duration  = numel(src.data) / o.quasiContinuousSamplingRate; % Ensure consistency
                    data = double(src.data')/1e6;% muV->V
                case "FIELDTRIP"
                    assert(o.duration<=src.hdr.nSamples/src.hdr.Fs,"Not enough samples in this FIELDTRIP struct");
                    o.recordingSamplingRate = round(src.fsample);
                    cfg = [];
                    cfg.channel= pv.channel;
                    cfg.latency = min(src.time{1})+[0 o.duration-1/src.hdr.Fs];
                    src = ft_selectdata(cfg,src); % Select channel and time window
                    data  = cat(2,src.trial{:})'/1e6;  % Make continuous data in V.
                case "RAW"
                    % A vector of doubles sampled at
                    % o.recordingSamplingRate
                    [nrSmples,nrChnls] = size(src); %#ok<ASGLU>
                    assert(nrSmples==o.recordingSamplingRate*o.duration,"The number of samples in this array (%d) does not match o.recordingSamplingRate*o.duration (%d)",nrSmples,o.recordingSamplingRate*o.duration);
                    data = src;
            end
            data = mean(data,2,"omitmissing"); % Average over channels
            % Fill missing  with linear interpolation (to avoid having nans)
            data = fillmissing(data,"linear");
            % Upsample to quasi continuous sampling rate to match the simulated data
            o.lfp = resample(data,o.quasiContinuousSamplingRate, o.recordingSamplingRate);
            o.lfp = o.lfp -mean(o.lfp,1,"omitmissing");
            o.fromData = true;
        end

        function sim(o)
            % Shortcut for the typical steps in a simulation.
            o.rng               % Reset rng if .reproducible
            o.reset;            % Start fresh
            o.simLfp;           % 1. Generate LFP/EEG .
            o.simSpikes;        % 2. Generate spikes
            o.simTCS;           % 3. Generate the TCS voltage
            o.simRecording;     % 4. Simulate recording
        end

        function simLfp(o)
            % simLfp - Simulate the Local Field Potential
            % LFPs can be single sinusoids (lfpPeakWidth =0) with
            % o.lfpAmplitude or be the sum of a distribution of lfps with an amplitude
            % that varies (as a Gaussian distribution with standard deviation lfpPeakWidth)
            % around o.lfpFrequency.
            %  PARMS used:
            %    lfpFrequency
            %    lfpAmplitude
            %    lfpPhaseOffset
            %    lfpPeakWidth
            %  COMPUTES
            %   lfp
            % If EEGLAB or Fieldtrip were used to set .lfp, skip the simulation.
            if o.fromData; return;end

            thisLfp = zeros(size(o.tSimulate));
            if all(o.lfpAmplitude ==0)
                o.lfp= thisLfp;
                return
            end

            nrOscillations= max([numel(o.lfpAmplitude) numel(o.lfpPeakWidth) numel(o.lfpPhaseOffset) numel(o.lfpFrequency)]);
            if isscalar(o.lfpAmplitude);o.lfpAmplitude = o.lfpAmplitude*ones(1,nrOscillations);end
            if isscalar(o.lfpPeakWidth);o.lfpPeakWidth= o.lfpPeakWidth*ones(1,nrOscillations);end
            if isscalar(o.lfpPhaseOffset);o.lfpPhaseOffset = o.lfpPhaseOffset*ones(1,nrOscillations);end
            if isscalar(o.lfpFrequency);o.lfpFrequency= o.lfpFrequency*ones(1,nrOscillations);end
            if any([numel(o.lfpAmplitude) numel(o.lfpPeakWidth) ...
                    numel(o.lfpPhaseOffset) numel(o.lfpFrequency)] ~= nrOscillations)
                error('All LFP parameter vectors must have the same length or be scalar.')
            end

            for i=1:nrOscillations
                if o.lfpPeakWidth(i) ==0
                    % Single oscillation .
                    thisLfp = thisLfp +  o.lfpAmplitude(i)*sin(2*pi*o.tSimulate*o.lfpFrequency(i)+o.lfpPhaseOffset(i));
                else
                    % Oscillation with frequency spread
                    % Construct Frequency Vector (0 to Nyquist)
                    if rem(o.nrTimePoints,2)==0
                        % Even number of samples. Nyquist is the highest one
                        frequencies = ([0 1:o.nrTimePoints/2 (o.nrTimePoints/2-1):-1:1]*o.quasiContinuousSamplingRate/o.nrTimePoints)';
                    else
                        % Odd number of samples. Nyquist is not in the set.
                        frequencies = ([0 1:(o.nrTimePoints-1)/2 (o.nrTimePoints-1)/2:-1:1]*o.quasiContinuousSamplingRate/o.nrTimePoints)';
                    end
                    [~,highestFrequencyIx] =max(frequencies);
                    % Select the lower half
                    keep = 1:highestFrequencyIx;
                    frequencies = frequencies(keep);

                    amplitude = o.lfpAmplitude(i) * exp(-0.5*((frequencies - o.lfpFrequency(i))/o.lfpPeakWidth(i)).^2);
                    phase = exp(1i * (2*pi*rand(size(frequencies)) - pi));
                    spectrum_half = amplitude .* phase;
                    spectrum_half(1) = 0; % Eliminate DC offset
                    if mod(o.nrTimePoints,2) == 0
                        spectrum_half(end) = real(spectrum_half(end));
                    end
                    % create the full two-sided spectrum
                    if mod(o.nrTimePoints,2) == 0
                        spectrum = [spectrum_half; conj(fliplr(spectrum_half(2:end-1)))];
                    else
                        spectrum = [spectrum_half; conj(fliplr(spectrum_half(2:end)))];
                    end
                    thisLfp = thisLfp + ifft(o.nrTimePoints/2*sqrt(frequencies(3)-frequencies(2))*spectrum, 'symmetric'); % 'symmetric' ensures real output
                end
            end
            o.lfp  = thisLfp-mean(thisLfp);
        end


        function simSpikes(o,pv)
            % simSpikes  - Generate a voltage trace representing the extracellular potential of spikes.
            %
            % Object properties used:
            % spikePeak
            % spikeRate
            % refractoryPeriod
            % poisson
            % preferredPhase - [rad]  (lock spikes to a certain phase)
            % synchrony  - Von Mises Kappa; large kappa means strong preference for the
            %           preferred phase.
            % fixNrSpikes - Set this to true to always generate the same number of spikes.
            %               Note that requiring a fixed number interacts with the
            %               goal to generate spikes at a certain phase but also a
            %               minimum refractoryPeriod.
            %
            % OPTIONAL Input parameters
            % kappa - Spike shape parameter representing the initial rise [5]
            % tau   - spike shape parameter representing the decay [5e-4]
            % graph  - generate a graph [false]
            % showTime  - limit the graph to the first n seconds of the stimulation [2]
            %
            % COMPUTES:
            % vSpike  = The voltage trace generated by the spikes only, sampled at the quasi-continuous rate [nrTimePoints 1]
            % ixSpike = The index into spikeV that corresponds to simulated spikes.
            % vSpikeNneighbor = The spiking voltage observed in a neighboring electtrode (same noise, same lfp, but a  different random draw of spikes)
            arguments
                o (1,1) artSim
                pv.kappa (1,1) double =5
                pv.tau (1,1) double =5e-5
                pv.graph (1,1) logical= false
                pv.showTime (1,1) double  =2
            end

            if o.spikePeak ==0
                % No spikes
                o.ixSpike =[];
                o.vSpike = zeros(o.nrTimePoints,1);
                o.vSpikeNeighbor=zeros(o.nrTimePoints,1);
                spikeTime = []; upDownSpike=[];
            else
                % Simulate one spike waveform
                spike= @(t,k,tau) (-1./k*tau*factorial(k-1))*(t./tau).^k.*exp(-t./tau);
                spikeTime = (0:1/o.quasiContinuousSamplingRate:o.refractoryPeriod-1/o.quasiContinuousSamplingRate)';
                bump  = spike(spikeTime,pv.kappa,pv.tau);
                upDownSpike = [0;diff(bump)];
                upDownSpike = o.spikePeak*upDownSpike./max(abs(upDownSpike)); %Scale to spikePeak voltage
                % Generat spikes twice. Once for the channel of interest, once for the
                % "neighboring" channel (different spike times, same lfp)
                for i=1:2
                    % Generate a spike train
                    nrSpk= round(o.duration*o.spikeRate);
                    if o.fixNrSpikes && nrSpk > floor(o.duration/o.refractoryPeriod)+1
                        error('fixNrSpikes cannot satisfy the requested rate and refractory period.')
                    end
                    if o.poisson
                        %Poisson
                        if o.synchrony >0
                            lambda = vonmises(o.phaseLfp,o.preferredPhase,o.synchrony);
                        else
                            lambda= ones(o.nrTimePoints,1);
                        end
                        lambda =lambda/sum(lambda)*(o.spikeRate*o.duration);
                        ix = find(poissrnd(lambda)>0);
                    else
                        % Generate nrSpikes to match rate using random spiking
                        ix = sort(randi(o.nrTimePoints,[nrSpk 1]));
                    end

                    % Avoid two spikes close together
                    out = [inf ;diff(ix)]/o.quasiContinuousSamplingRate < o.refractoryPeriod;
                    ix(out)=[];

                    % Remove within refractoryPeriod and add spikes until we get to the target
                    % number. Note that this can remove any phase concentration (synchrony>0) that
                    % may have been created above (because spikes at the preferred phase are
                    % more likely to be within refractoryPeriod from each other. So for large
                    % kappa and high firing rates, this will produce weird results
                    while o.fixNrSpikes
                        % Remove
                        out = [inf ;diff(ix)]/o.quasiContinuousSamplingRate< o.refractoryPeriod;
                        ix(out)=[];

                        % Check whether we have too few or too many spikes and remove/add to
                        % match the desired mean spike rate
                        nrToGenerate = nrSpk-numel(ix);
                        if nrToGenerate<0
                            out= randperm(numel(ix));
                            ix(out(1:abs(nrToGenerate)))=[];
                        elseif nrToGenerate>0
                            if o.poisson
                                extra =  find(poissrnd(lambda)>0);
                            else
                                extra = sort(randi(o.nrTimePoints,[nrToGenerate 1]));
                            end
                            take = randperm(numel(extra));
                            ix =[ix;extra(take(1:nrToGenerate))];
                            ix = sort(ix);
                        else %nrToGenerate==0
                            break; % Step out of while
                        end
                    end

                    % Generate the voltage trace by convolution
                    isSpike= zeros(o.nrTimePoints,1);
                    isSpike(ix)=1+o.spikePeakVariation*randn([numel(ix),1]);

                    V = conv(isSpike,upDownSpike,'full');
                    V(end-numel(spikeTime)+2:end)=[];
                    if i==1
                        % The first one is the neighbor; for these ixSpike is not used and
                        % not returned.
                        o.vSpikeNeighbor = V;
                    else
                        o.vSpike = V;
                        o.ixSpike = ix;
                    end
                end
            end

            %% Show results
            if pv.graph
                clf;
                subplot(2,1,1)
                plot(spikeTime/1e-3,upDownSpike/1e-6)
                xlabel 'Time (ms)'
                ylabel 'Voltage (\muV)'
                title 'Single Spike Waveform'

                subplot(2,1,2);
                plot(o.tSimulate,o.vNeural/1e-6)
                hold on
                plot(o.tSimulate(o.ixSpike),zeros(size(o.ixSpike)),'*');
                xlim([0 pv.showTime]);
                xlabel 'Time (s)'
                ylabel 'Voltage (\muV)'
                title 'Trace'
            end
        end

        function simTCS(o,pv)
            % simTCS - compute the stimulation voltage
            %
            % Object properties used:
            % tacsAmplitude  - Amplitude of the sinusoid [A]
            % tacsFrequency   - Frequency of stimulation [Hz]
            % tacsPhaseOffset     - Phase at time 0. [rad]
            % tacsShape    - sine,sawtooth, square, am
            % tacsCarrierFrequency  - Carrier frequency of the amplitude-modulated tACS
            % tdcsMean   -  Mean current [A]
            % tcsRamp   - Duration of a ramp [s]
            %
            %
            % resistance - Effective resistance of the path between stimulation
            %               electrode and the recording site.
            % stimulatorSamplingRate - Sampling rate of the stimulation device [Hz]
            % currentResolution - Resolution of the stimulator [A]
            %
            % OPTIONAL INPUT PARMS
            % graph     = Toggle to show a grpah
            % showTime - Show the first seconds of stimulation in the graph [2]
            %
            % COMPUTES
            % vTcsActual =  Voltage as generated by the stimulator (sampled at .stimulatorSamplingRate)

            arguments
                o (1,1) artSim
                pv.graph (1,1) logical =false
                pv.showTime (1,1) double = 2
            end

            assert(o.stimulatorSamplingRate>max(o.tacsFrequency,o.tacsCarrierFrequency)*2,"The stimulator sampling rate (%f) is too low for this carrier frequency (%f).",o.stimulatorSamplingRate,o.tacsCarrierFrequency)
            % Determine the output
            switch upper(o.tacsShape)
                case 'SINE'
                    virtualOutputCurrent = o.tacsAmplitude*sin(2*pi*o.tacsFrequency*o.tTcs+o.tacsPhaseOffset);
                case 'SAWTOOTH'
                    % Sawtooth - with sine-phase (i.e. 0 at 0).
                    virtualOutputCurrent = o.tacsAmplitude*sawtooth(2*pi*o.tacsFrequency*o.tTcs+o.tacsPhaseOffset+pi);
                case 'SQUARE'
                    virtualOutputCurrent = o.tacsAmplitude*(2*(sin(2*pi*o.tacsFrequency*o.tTcs+o.tacsPhaseOffset)>0)-1);
                case 'NONE'
                    virtualOutputCurrent = zeros(size(o.tTcs));
                case 'AM'
                    % Amplitude modulated sinusoidal tACS
                    virtualOutputCurrent = o.tacsAmplitude*(0.5*sin(2*pi*o.tacsFrequency*o.tTcs+o.tacsPhaseOffset)+0.5).*sin(2*pi*o.tacsCarrierFrequency*o.tTcs);
                otherwise
                    error('Unknown tACS shape %s',o.tacsShape);
            end

            if o.tdcsMean~=0
                virtualOutputCurrent = virtualOutputCurrent + o.tdcsMean*ones(numel(o.tTcs),1);
            end
            if o.tcsRamp >0
                nrRampSamples = min(floor(numel(o.tTcs)/2),sum(o.tTcs <= o.tcsRamp));
                scale = ones(size(o.tTcs));
                if nrRampSamples > 0
                    scale(1:nrRampSamples) = linspace(0,1,nrRampSamples)';
                    scale(end-nrRampSamples+1:end) = linspace(1,0,nrRampSamples)';
                end
            else
                scale =1;
            end
            % Model the discretization of the DA conversion in the device
            I = o.currentResolution*round(scale.*virtualOutputCurrent./o.currentResolution);

            % The stimulator applies this current to the scalp and that ultimately results
            % in a voltage difference between the intracranial electrode near the neuron of
            % interest, and a reference electrode. For simplicity we assume that this involves
            % no nonlinear changes (i.e. everything is purely ohmic).
            % Store the voltage.
            V = I.*o.resistance;

            % The stimulator updates its value with stimulatorSamplingRate; values are
            % constant in between. Create a v that reflects this discrete
            % nature of the stimulation
            o.vTcsActual= interp1(o.tTcs,V,o.tSimulate,"previous","extrap");

            %% Simulate impedance artifacts
            if isa(o.zFun,'function_handle')
                % User-supplied function, pass the object
                z = o.zFun(o);
            elseif o.zAmplitude >0
                % Simulate a periodic breathing or heartbeat artifact.
                if o.zFrequency <= 0
                    z = ones(size(o.tSimulate));
                else
                    if o.tacsFrequency>0 && o.tacsFrequency/o.zFrequency == round(o.tacsFrequency/o.zFrequency)
                        fprintf(2,'Heartbeat is effectively locked to tACS - this is not a good simulation of a heartbeat... (use a frequency that is not a multiple of the tACS frequency)\n');
                    end
                    % Simulate a periodic resistance change with optional long-term variability.
                    frequency = o.zFrequency*(1+o.zVariability*sin(2*pi*(1/3600)*o.tSimulate));
                    z = sin(2*pi*frequency.*o.tSimulate);
                    if o.zDuration >0
                        % Create bumps with a given FWHM by raising the sine to a power of n.
                        n = (-log(2)/log(cos(pi*o.zFrequency*o.zDuration)));
                        z(z<0) = 0;
                        z= z.^n;
                        z = z./max(z);
                    end
                    z = 1+ o.zAmplitude.*z;
                end
            else
                z= 1;
            end

            o.vTcsActual = o.vTcsActual.*z;



            %% Show a summary graph
            if pv.graph
                clf;
                subplot(3,1,1)
                plot(o.tTcs,I*1000)
                xlabel ('Time (s)')
                ylabel 'Current (mA)'
                title('Stimulator Output Current');
                subplot(3,1,2)
                plot(o.tTcs,I*1000)
                xlim([0 pv.showTime]);
                xlabel ('Time (s)')
                ylabel 'Current (mA)'
                title('Zoomed view');

                subplot(3,1,3)
                [ft,frequency]= fftReal(o.vTcsActual(:),o.quasiContinuousSamplingRate);
                amplitude = abs(ft)/(numel(o.vTcsActual)/2);
                amplitude(1)= 0; % Remove DC
                [maxAmp,maxIx] = max(amplitude);
                plot(frequency,amplitude./maxAmp);
                xlabel 'Frequency (Hz)'
                ylabel 'Relative Amplitude'
                title(sprintf('Amplitude Spectrum - Max (%3.2g mV @ %3.3g Hz)',maxAmp/1e-3,frequency(maxIx)))
            end
        end

        function  simRecording(o,pv)
            %simRecording - Simulate the output of the amplifier that records
            % the combination of a stimulation artifact and neural signals.
            %
            % OPTIONAL EXTRA PARM/VALUE PAIRS
            % graph - Show a graph  [false]
            % showTime - limit graph to the first seconds [2]
            %
            % COMPUTES
            %  vRecord - The Voltage as recorded by the device. The first column represents
            %       the voltage during sham (including  additive noise ). The second column
            %       represents the same neural signal and noise,now with the stimulation artifacts (including multiplicative noise).
            % vRecordNeighbor - Same as vRecord, but for a neihboring
            %                       electrode (Same lfp , different
            %                       spikes).
            %
            %  vRecordTcs- the applied stimulation voltage, as recorded by the
            %               recording device (and thus reflecting
            %               digitization of the stimulator).
            arguments
                o (1,1) artSim
                pv.graph (1,1) logical =false
                pv.showTime (1,1) double = 2
            end
            
            %% Record the TCS signal coming from the stimulator.
            % Assume a linear amplifier ; only downsampling (includes
            % anti-aliasing filter in decimate).
            o.vRecordTcs = decimate(o.vTcsActual,o.quasiContinuousSamplingRate/o.recordingSamplingRate);

            % Create 4 recorded signals; one for the recording electrode,
            % one for a neighboring electrode, both with and without tCS.
            for i=1:2  % neighbor (1), electrode(2)
                if i==1
                    simV = o.vNeuralNeighbor;
                else
                    simV = o.vNeural;
                end
                % Additive noise
                if o.pinkNoise
                    additiveNoise= o.additive*coloredNoise(nrTimePoints=o.nrTimePoints,alpha=1);
                else
                    additiveNoise = o.additive*randn(o.nrTimePoints,1);
                end
                simV = simV + additiveNoise;                
                simV = simV.*(1+o.multiplicative.*randn(o.nrTimePoints,1));

                for j=1:2  % 1 -> without tacs (.vTruth), 2-> with tACS
                    if j==2
                        simV = simV +o.vTcsActual.*(1+o.multiplicative.*randn(o.nrTimePoints,1)); % Add the actual tcs signal
                    end
                    % Simulate the preamp/amp/post-amp filter chain
                    recV = adc(o,simV);
                    if i==1
                        o.vRecordNeighbor(:,j) = recV;
                    else
                        o.vRecord(:,j) = recV;
                    end
                end
            end

            %% Summarize the results
            if pv.graph
                range = prctile(o.vRecord,[0.5 99.5])/1e-3;
                bit = o.vRecord/o.adcVoltsPerBit;
                headroom = (2.^o.adcBitDepth-1) - max(abs(bit),[],1);

                for i=1:2
                    subplot(2,2,(i-1)*2+1)
                    plot(o.tRecord,o.vRecord(:,i)/1e-3)
                    xlabel ('Time (s)')
                    ylabel 'Voltage (mV)'
                    title(sprintf('V99%%: [%3.2f %3.2f] Headroom %d ADU',range(1,i),range(2,i)),headroom(i));
                    xlim([0 pv.showTime])

                    subplot(2,2,(i-1)*2+2)
                    [ft,frequency]= fftReal(o.vRecord(:,i),o.recordingSamplingRate);
                    amplitude = abs(ft)/(numel(o.vRecord(:,1))/2);
                    [maxAmp,maxIx] = max(amplitude);
                    plot(frequency,amplitude/1e-3);
                    xlabel 'Frequency (Hz)'
                    ylabel 'Amplitude (mV)'
                    set(gca,'XScale','Log')
                    title(sprintf('Spectrum - Max (%3.2g mV @ %3.3g Hz)',maxAmp/1e-3,frequency(maxIx)))


                    if i==1
                        ax = axes('Position',[.8 .7 .1 .1]);
                        % Show inset with amplification profile plus saturation range.
                        x = min(o.vTcsActual(:)):1e-6:max(o.vTcsActual(:));
                        amplifierInputRange = o.amplifierLinearRange/o.amplifierGain;
                        [y,satPoint]= artSim.amplify(x,amplifierInputRange,o.adcNonlinearity,o.amplifierGain,satLevel=99);
                        plot(ax,x/1e-3,y/1e-3);
                        hold on
                        plot(satPoint*[1e3 1e3],ylim)
                        plot(-satPoint*[1e3 1e3],ylim)
                        xlabel 'mV'
                        ylabel 'mV'
                    end
                end
            end

        end

        function [v,saturationPoint] = adc(o,v,pv)
            % Take a quasi continuous signal and simulate the preamp, amp,
            % post-amp filtering stages
            arguments
                o (1,1) artSim
                v (:,1) double
                pv.satLevel (1,1) double = 99;
                pv.isTimeSeries (1,1) logical = true;
            end

            if o.highPass >0 && pv.isTimeSeries
                highPassFiltered  = highpass(v,o.highPass,o.quasiContinuousSamplingRate);
            else
                highPassFiltered  = v;
            end
            % Analog preamplifier gain
            amplifierInputRange = o.amplifierLinearRange/o.amplifierGain;
            [amplifiedV,saturationPoint] = artSim.amplify(highPassFiltered, ...
                amplifierInputRange,o.adcNonlinearity,o.amplifierGain,satLevel=pv.satLevel);
            % Analog anti-alias filter before sampling.
            if o.lowPass >0 && pv.isTimeSeries
                amplifiedV = lowpass(amplifiedV,o.lowPass, ...
                    o.quasiContinuousSamplingRate);
            end
            % Sampling
            if pv.isTimeSeries
                if o.lowPass >0
                    sampledV = downsample(amplifiedV, ...
                        o.quasiContinuousSamplingRate/o.recordingSamplingRate);
                else
                    % decimate supplies the anti-aliasing filter when lowPass is zero.
                    sampledV = decimate(amplifiedV, ...
                        o.quasiContinuousSamplingRate/o.recordingSamplingRate);
                end
            else
                sampledV = amplifiedV;
            end
            sampledV = min(max(sampledV,-o.adcRange),o.adcRange);
            % ADC quantization.
            v = o.adcVoltsPerBit *round(sampledV./o.adcVoltsPerBit);
            % Optional digital filter after ADC quantization.
            if o.digitalLowPass >0 && pv.isTimeSeries
                v = lowpass(v,o.digitalLowPass, o.recordingSamplingRate);
            end
        end

        function [results,spikes] = sortSpikes(o, vClean,pv)
            % [results,spikes] = sortSpikes(o, vClean,pv)
            % Simulate a typical spike sorting pipelinw.
            % First filter the data to the 'spikeBand', then use
            % UMS2K to detect spikes based on a 'threshold', and keep
            % either 'ALL' detected spikes, the cluster with 'MOSTSPIKES',
            % or the cluster that matches the waveforms that occur at least 'minNrSpikes' times
            % in the sham (unstimulated) trials ('MATCHSHAM')
            % INPUT
            % 'threshold' - Standard deviations of the signal [3.5]
            % 'sharedThreshold' - Set to true to determine the threshold by
            %                   pooling all data (tacs+ sham), set to false to use a separate
            %                   voltage threshold for tacs and sham trials.
            %                   [true]
            % 'sortMode' = MATCHSHAM, MOSTSPIKES,ALL
            % 'spikeBand' - Frequency band in which to detect spikes [300 7000] Hz.
            % 'filterOrder' - Order of the Butterworth Bandpass filter [3].
            % 'minNrSpikes' - A spike count above this in the Sham trial is
            %                   considered a neuron. (used only by
            %                   MATCHSHAM)
            % phaseFromRaw - Set to true to determine phase from the
            %                   recorded signal, false uses the applied
            %                   tACS
            % nrPhaseBins  - Number of bins for phase.
            % nrBootstraps - Bootstrap sets for significance estimate of
            %                   PLV
            %
            % Output
            % results - Struct array with results for the sham trials (1)
            % and for trials with tacs (2)
            % spikes - Raw output of UMS2K.
            %
            arguments
                o (1,1) artSim
                vClean (:,1) double
                pv.threshold (1,1) double =3.5
                pv.sharedThreshold (1,1) logical = true
                pv.sortMode (1,1) string {mustBeMember(pv.sortMode,["MATCHSHAM","MOSTSPIKES" "ALL" ])} = "MATCHSHAM"
                pv.spikeBand (1,2) double = [300 7000]
                pv.filterOrder (1,1) double = 3
                pv.minNrSpikes (1,1) double = 100
                pv.phaseFromRaw (1,1) logical = true
                pv.nrPhaseBins (1,1) double = 8
                pv.nrBootstraps (1,1) double  =0
            end

            % Filter the signal in the band that contains the spikes.
            % A Butterworth bandpass filter is a common way to do this:
            if all(isfinite(pv.spikeBand))
                [A,B] = butter(pv.filterOrder,pv.spikeBand/(o.recordingSamplingRate/2));
                filteredVoltage = filtfilt(A,B,[o.vTruth vClean ]);
            else
                filteredVoltage = [o.vTruth vClean ];
            end


            %% Do spike detection and spike sorting using UMS2K
            spikes = ss_default_params(o.recordingSamplingRate);
            spikes.params.window_size = 2+1000/o.recordingSamplingRate;
            spikes.params.cross_time = 1;
            spikes.params.thresh = pv.threshold;
            if pv.sharedThreshold
                % Determine one threshold for tACS and Sham trials together
                spikes = ss_detect({filteredVoltage(:,1),filteredVoltage(:,2)},spikes); % Trial  1= ground truth, Trial 2 = recovered signal
            else
                % Estimate separate thresholds
                % 1. ground truth
                spikes.params.thresh =  std(filteredVoltage(:,1))* -pv.threshold;
                spikes.params.detect_method ='manual';
                shamSpikes = ss_detect({filteredVoltage(:,1)},spikes);

                %Trial 2 = recovered signal
                spikes.params.thresh =  std(filteredVoltage(:,2))* -pv.threshold;
                spikes.params.detect_method ='manual';
                stimulatedSpikes = ss_detect({filteredVoltage(:,2)},spikes);

                % Combine the independently detected trials before alignment and
                % sorting.
                spikes = shamSpikes;
                nrShamSpikes = numel(shamSpikes.spiketimes);
                nrStimulatedSpikes = numel(stimulatedSpikes.spiketimes);
                spikes.spiketimes = [shamSpikes.spiketimes stimulatedSpikes.spiketimes];
                spikes.unwrapped_times = [shamSpikes.unwrapped_times stimulatedSpikes.unwrapped_times];
                spikes.waveforms = [shamSpikes.waveforms; stimulatedSpikes.waveforms];
                spikes.trials = [ones(1,nrShamSpikes)  2*ones(1,nrStimulatedSpikes)];
                spikes.info.detect.event_channel = ones(nrShamSpikes+nrStimulatedSpikes,1);
                spikes.info.detect.thresh = [shamSpikes.info.detect.thresh stimulatedSpikes.info.detect.thresh];
            end

            spikes = ss_align(spikes);
            spikes = ss_kmeans(spikes);
            spikes = ss_energy(spikes);
            spikes = ss_aggregate(spikes);
            isStimulated = spikes.trials ==2;

            switch upper(pv.sortMode)
                case 'MATCHSHAM'
                    % Do spike sorting; keep spikes in the stimulated trial that look like the spikes that
                    % occur at least minNrSpikes times in the unstimulated trials
                    [spikeClass,~,ic]= unique(spikes.assigns(~isStimulated));
                    nrSpikesPerClass = accumarray(ic,1);
                    stay  = nrSpikesPerClass>pv.minNrSpikes;
                    theNeuron = spikeClass(stay);
                case 'MOSTSPIKES'
                    % Keep only the unit with the most spikes
                    [spikeClass,~,ic]= unique(spikes.assigns);
                    nrSpikesPerClass = accumarray(ic,1);
                    [~,ix] = max(nrSpikesPerClass);
                    theNeuron = spikeClass(ix);
                case 'ALL'
                    % Keep all threshold crossing (i.e. not spike sorting)
                    theNeuron = unique(spikes.assigns);
                otherwise
                    error('Unknown spikeSortingMode %s',sortMode)
            end

            %% Post process
            % A found spike is a spike within +/- 1 ms from a simulated spike
            withinOneMs = @(x,y)(ismember(x,y) | ismember(x-1,y) | ismember(x+1,y));

            ixDetectedSpike = round(spikes.spiketimes*o.recordingSamplingRate);
            detectedSpikeMs = round(spikes.spiketimes*1000);
            spikeMs         = round(1000*o.ixSpike/o.quasiContinuousSamplingRate);

            % Store
            for i = 1:2  %Sham /TACS
                if i==1
                    stay = ~isStimulated;
                    label = 'SHAM';
                else
                    stay = isStimulated;
                    label ='tACS';
                end
                isAssignedSpike = stay & ismember(spikes.assigns,theNeuron);
                results(i).ixDetected = ixDetectedSpike(isAssignedSpike); %#ok<*AGROW>
                results(i).detectedMs = detectedSpikeMs(isAssignedSpike);
                results(i).detectedMs(results(i).detectedMs ==0) =1;% Used as an index below ; must be>0.

                results(i).WF         = spikes.waveforms(isAssignedSpike,:)';
                % Detected but not assigned as "theNeuron":
                results(i).unassignedWF   = spikes.waveforms(stay & ~ismember(spikes.assigns,theNeuron),:)';

                results(i).found      = withinOneMs(spikeMs,results(i).detectedMs);
                results(i).invented   = ~withinOneMs(results(i).detectedMs,spikeMs);
                results(i).inventedWF = results(i).WF(:,results(i).invented);

                results(i).ixMissed  =  round(o.ixSpike(~results(i).found)*o.recordingSamplingRate/o.quasiContinuousSamplingRate);
                results(i).missedMs  =  round(1000*o.ixSpike(~results(i).found)/o.quasiContinuousSamplingRate);
                results(i).missedMs(results(i).missedMs  ==0) =1;% Used as an index below ; must be>0.
                nrMissed = numel(results(i).missedMs);
                if nrMissed>0
                    % Extract waveforms of missed spikes using the same
                    % approach as UMS2K
                    window_samples  = round( spikes.params.Fs * spikes.params.window_size / 1000);
                    samples_before  = round( spikes.params.Fs * spikes.params.cross_time /1000);
                    samples_after   = round( spikes.params.Fs * spikes.params.max_jitter / 1000)+ window_samples - (1+samples_before);
                    nrSamplesPerWF = samples_after+samples_before+1;
                    ix = results(i).ixMissed +  (-samples_before:samples_after);
                    out = ix<1 | ix > size(filteredVoltage,1);
                    ix(out) = 1;
                    results(i).missedWF =reshape(filteredVoltage(ix',i),nrSamplesPerWF,[]);

                    % Hack to use UMS2K code for alignment of the
                    % missed spikes (to match detected spikes).
                    tmpSpikes =spikes;
                    tmpSpikes.info = rmfield(tmpSpikes.info,'align');
                    tmpSpikes.waveforms = results(i).missedWF';
                    tmpSpikes.spiketimes = results(i).missedMs'/1000; % seconds
                    tmpSpikes.info.detect.event_channel = ones(nrMissed,1);
                    tmpSpikes= ss_align(tmpSpikes);
                    results(i).missedWF=tmpSpikes.waveforms';
                else
                    results(i).missedWF = [];
                end
                results(i).voltage    = filteredVoltage(:,i);
                results(i).label = label;
            end



            %% Analyze waveforms with respect to tACS Phase.
            % Downsample the voltage to 1 kHz to speed up filtering
            assert(o.tacsFrequency<500,"Adjust filtering to analyze high tACS frequencies");
            if pv.phaseFromRaw
                % Estimate the phase from the rawVoltage on the neighboring electrode
                % to avoid spike bleed.
                voltageInBand = resample(o.vRecordNeighbor,1000,o.recordingSamplingRate);
            else
                % Estimate phase from the tACS signal.
                voltageInBand = resample(o.vRecordTcs,1000,o.recordingSamplingRate);
                voltageInBand = repmat(voltageInBand,[1 2]);
            end
            if o.tacsFrequency>0
                [B,A] = butter(3,(o.tacsFrequency+[-1 +1])/(1000/2),'bandpass');
                voltageInBand = filtfilt(B,A,voltageInBand );
            end
            estimatedPhase = artSim.phase(voltageInBand);
            estimatedPhase = mod(estimatedPhase+2*pi,2*pi); % Warp to [0 2pi];
            % Determine PPC/PLV
            for i=1:2
                results(i).estimatedPhase = estimatedPhase(:,i);
                [results(i).ppc,results(i).pPpc] = ppc(results(i).detectedMs,results(i).estimatedPhase,pv.nrBootstraps);
                [results(i).ppc2,results(i).pPpc2] = ppc(results(i).detectedMs,2*results(i).estimatedPhase,pv.nrBootstraps);
            end

            % Bin waveforms relative to the "lfp" or tACS phase
            binWidth = 2*pi/pv.nrPhaseBins;
            for i=1:2
                results(i).meanWF = nan(size(spikes.waveforms,2),pv.nrPhaseBins);
                results(i).sdWF = nan(size(spikes.waveforms,2),pv.nrPhaseBins);
                results(i).meanMissedWF = nan(size(spikes.waveforms,2),pv.nrPhaseBins);
                results(i).sdMissedWF = nan(size(spikes.waveforms,2),pv.nrPhaseBins);


                edges= linspace(0,2*pi,pv.nrPhaseBins+1);
                results(i).binCenters = edges(1:end-1)';
                phaseBins =  discretize(mod(results(i).estimatedPhase+0.5*binWidth, 2*pi), edges);
                results(i).phaseCount = accumarray(phaseBins,ones(size(phaseBins)),[],@sum);
                % For the found spikes.
                if ~isempty(results(i).WF)
                    phaseBinPerSpike = phaseBins(results(i).detectedMs); % phase is in 1kHz so detectedMs can be used as index
                    results(i).spikeCount = accumarray(phaseBinPerSpike,ones(size(phaseBinPerSpike)),[],@sum);
                    for j=1:pv.nrPhaseBins
                        stay = phaseBinPerSpike ==j;
                        if any(stay)
                            results(i).meanWF(:,j) = mean(results(i).WF(:,stay),2);
                            results(i).sdWF(:,j) = std(results(i).WF(:,stay),1,2);
                        end
                    end
                else
                    results(i).spikeCount  =zeros(pv.nrPhaseBins,1);
                end
                % Same for the missed spikes
                if ~isempty(results(i).missedWF)
                    phaseBinPerSpike = phaseBins(results(i).missedMs); % phase is in 1kHz so detectedMs can be used as index
                    results(i).missedSpikeCount = accumarray(phaseBinPerSpike,ones(size(phaseBinPerSpike)),[],@sum);
                    for j=1:pv.nrPhaseBins
                        stay = phaseBinPerSpike==j;
                        if any(stay)
                            results(i).meanMissedWF(:,j) = mean(results(i).missedWF(:,stay),2);
                            results(i).sdMissedWF(:,j) = std(results(i).missedWF(:,stay),1,2);
                        end
                    end
                else
                    results(i).missedSpikeCount  =zeros(pv.nrPhaseBins,1);
                end
                % accumarry does not count phase bins that never occured.
                % Set to 0.
                nrMissing = pv.nrPhaseBins-numel(results(i).phaseCount);
                results(i).phaseCount = [results(i).phaseCount; zeros(nrMissing,1)];
                results(i).spikeCount = [results(i).spikeCount; zeros(nrMissing,1)];
                results(i).missedSpikeCount = [results(i).missedSpikeCount; zeros(nrMissing,1)];
            end
        end
    end

    methods (Static)
        function v = phase(x)
            % Determines the (sine) phase of a signal using the Hilbert
            % transform - add 0.5 pi to get sin phase
            v = angle(hilbert(x))+0.5*pi;
        end

        function [y,saturationPoint] = amplify(x,range,n,gain,pv)
            % [y,saturationPoint] = amplify(x,range,n,gain,pv)
            % Simulated signal amplification with a specified gain factor,
            % a linear range and adjustable saturation.
            %
            % For large n (e.g., 100), the function is piecewise linear, with linear amplification
            % up to +/- range and saturation at the range for larger/smaller values.
            %
            % For small n (e..g, 2) the function is sigmoidal; the saturation pont (at
            % any percentage of the input) can be computed by providing 'satLevel'
            %
            % INPUT
            %  x = signal
            % range = linear range of the amplifier.
            % n = nonlinearity of the amplifier.
            % gain -  Gain factor. Defaults to 1.
            % Parameter/Value pairs
            % satLevel = Compute x-value for which amplification is this percentage of
            %               the input. Default 99 corresponds to "linear within 1%".
            % OUTPUT
            % y = amplified signal
            %
            % EXAMPLE
            %
            %             x= -10:0.1:10; % Input
            %             range = 6.8;  % Nominally linear range
            %             n = 20;       % Nonlinearity
            %             satLevel = 99; % Estimate nonlinearity at the 99% level.
            %
            %             [y,satPoint] = amplify(x,range,n,1,satLevel=satLevel);
            %             clf
            %             plot(x,y);
            %             hold on
            %             plot(satPoint*[1 1],ylim)
            %             xlabel 'Input'
            %             ylabel 'Output'
            %             title(sprintf('Range %3.1f Nonlinearity %d is %3.1f %% linear up to %3.2f',range,n,satLevel,satPoint))
            %
            arguments
                x           double
                range       double
                n           double {mustBeInteger,mustBeGreaterThan(n,1)}
                gain        double  = 1
                pv.satLevel double  = 99
            end

            sX =sign(x);
            aX = abs(x);
            y= gain.*sX.*(aX./(1+(aX./range).^n).^(1/n));
            saturationPoint = ((100/pv.satLevel)^n-1)^(1/n)*range;
        end

    end
end
