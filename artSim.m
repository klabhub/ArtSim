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
        % Set it to a nonnegative integer to use that as a seed.

        duration                  = 120;       % [s] Duration of the recording to simulate.

        recordingSamplingRate     = 30e3;     % [Hz] Sampling rate of the recording device
        stimulatorSamplingRate    = 2e3;      % [Hz] Sampling rate of the tACS stimulator
        currentResolution         = 1e-6;     % [A] Current resolution of the stimulator.
        highPass                  = 0;        % [Hz] Highpass cutoff of the pre-amplification hardware filter.

        % tACS
        tacsAmplitude             = 0;        % [A] Current amplitude
        tacsPhaseOffset           = 0;        % [Rad] Current phase
        tacsFrequency             = 0;        % [Hz] Current frequency
        tacsShape                 = 'sine';   % ['sine'] {'sine','sawtooth','square'}
        % tDCS
        tdcsMean                  = 0;       % [A] Current
        % Ramp (applied to both tACS and tDCS)
        tcsRamp                  = 0;       % [s] Duration of the on and offset ramp.

        resistance                = 1;        % [Ohm] Effective resistance at the recording site.

        % Simulating periodic resistance changes as caused by physiological changes such as heartbeats or respiration
        % The default implementation uses a periodic sinusoidal or
        % pulse-like shape, but the user can specify any zFun to simulate
        % more complex shapes (see zHeartbeat.m for an example)
        zFrequency        = 0;        % [Hz]
        zAmplitude        = 0;        % [fraction of the resistance]
        zDuration         = 0.25      % [s] - Full width at half max.
        zVariability      = 0;        % [fraction of heartbeat frequency per *hour*]
        zFun             = []         % Function that returns the impedance at each timepoint.

        % Spikes
        spikePeak                 = 0;        % [V] Peak of the action potential
        spikePeakVariation        = 0;        % Spikes are scaled by 1+randn*spikePeakVariation.
        spikeRate                 = 0;        % [Hz] Firing rate
        poisson                   = true;     % Set to true to generate spike times with a Poission distribution .
        % When false, spike times are drawn from a flat distribution.
        fixNrSpikes               = false;    % Fix the number of spikes to match the rate. (see simSpikes)
        refractoryPeriod          = 1.5e-3;   % [s] Refractory period
        synchrony                 = 0;        % Von Mises Kappa:
        % 0 means Poisson random spike times,
        % Typical PLV <0.05 which corresponds to a synchrony of 0.1
        % A high PLV of 0.3 is reached for synchrony=1;
        lfpAmplitude              = 0;        % [V] amplitude of the intrinsic LFP signal.
        lfpPhaseOffset            = 0;        % [rad]
        lfpFrequency              = 0;        % [Hz] Frequency of an intrinsic oscillation
        lfpPeakWidth              = 0;        % [Hz] Standard deviation of the ampltiude distribution for a collection of LFPs with random phases.
        preferredPhase            = 0;        % [Rad] Phase at which spikes occur preferentially (if syncrony>0)

        additive                  = 0;        % [V]  Standard deviation of additive noise
        pinkNoise                 = true;     % Set to true to generate pink noise (white otherwise)
        multiplicative            = 0;        % Standard deviation of the multiplicative noise

        %% Properties of the ADC converter.
        adcLinearRange            = 6.4e-3;     % [V] Range over which the amplifier is linear
        adcNonlinearity           = 100;       % The saturation profile of the amplifier near the edge
        %                                       of the linear range. 2 = very sigmoidal, 100 = linear inside the range,
        %                                       saturation outside (i.e. piecewise linear). See amplify
        adcRange                   = 10e-3;     % The range of the ADC; bits are distributed evenly across this +/- range.
        adcBitDepth                = 16;       % Number of bits in the ADC
        adcGain                    = 1;        % Amplifier gain
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
        nrSamples;  % The number of simulated recorded samples
        nrTimePoints; % The number of simulated time points
        phaseLfp;   % The phase of the LFP at each time point.
        vNeural;    % The neural voltage trace (spikes + lfp)
        vNeuralNeighbor   % The neural voltage trace at a neighboring electrode.
        vTruth;     % The true voltage (as would have been recorded without the artifact)
        vContaminated; % The recorded voltage (with the artifact).
        quasiContinuousSamplingRate;  % "Sampling" rate of the underlying simulation
        phaseTacs; % The phase of the applied/intended tACS signal
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
        function reset(o)
            % Reset the object.
            o.ixSpike = [];
            o.vSpike = zeros(o.nrTimePoints,1);
            o.vSpikeNeighbor = zeros(o.nrTimePoints,1);
            o.lfp = zeros(o.nrTimePoints,1);
            o.tRecord = (0:1/o.recordingSamplingRate:(o.duration-1/o.recordingSamplingRate))';
            o.tSimulate = (0:1/o.quasiContinuousSamplingRate:(o.duration-1/o.quasiContinuousSamplingRate))';
            o.vRecord = zeros(o.nrSamples,2);
            o.vRecordNeighbor = zeros(o.nrSamples,2);
            o.vRecordTcs = zeros(o.nrSamples,1);
        end
        function o = artSim(prms)
            %artSim Construct an artSim instance
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
            %   lfpPhase
            if o.lfpAmplitude ==0;return;end
            thisLfp = zeros(size(o.tSimulate));
            nrOscillations= max([numel(o.lfpAmplitude) numel(o.lfpPeakWidth) numel(o.lfpPhaseOffset) numel(o.lfpFrequency)]);
            if isscalar(o.lfpAmplitude);o.lfpAmplitude = o.lfpAmplitude*ones(1,nrOscillations);end
            if isscalar(o.lfpPeakWidth);o.lfpPeakWidth= o.lfpPeakWidth*ones(1,nrOscillations);end
            if isscalar(o.lfpPhaseOffset);o.lfpPhaseOffset = o.lfpPhaseOffset*ones(1,nrOscillations);end
            if isscalar(o.lfpFrequency);o.lfpFrequency= o.lfpFrequency*ones(1,nrOscillations);end

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


                    amplitude = o.lfpAmplitude(i) * exp(-((frequencies - o.lfpFrequency(i))/o.lfpPeakWidth(i)).^2);
                    phase = exp(1i * (2*pi*rand(size(frequencies)) - pi));
                    spectrum_half = amplitude .* phase;
                    % Enforce Real Signal Constraints (DC must be real, Nyquist real)
                    spectrum_half(1) = 0; % Eliminate DC offset
                    if mod(o.nrTimePoints,2) == 0
                        spectrum_half(end) = real(spectrum_half(end));
                    end
                    % G. Create Full Two-Sided Spectrum (Conjugate Symmetry)
                    if mod(o.nrTimePoints,2) == 0
                        spectrum = [spectrum_half; conj(fliplr(spectrum_half(2:end-1)))];
                    else
                        spectrum = [spectrum_half; conj(fliplr(spectrum_half(2:end)))];
                    end
                    thisLfp = thisLfp + ifft(o.nrTimePoints/2*sqrt(frequencies(3)-frequencies(2))*spectrum, 'symmetric'); % 'symmetric' ensures real output
                end
            end
            o.lfp  = thisLfp;
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
            % graph  - generate a graph [true]
            % showTime  - limit the graph to the first n seconds of the stimulation [2]
            %
            % COMPUTES:
            % vSpikes  = The voltage trace generated by the spikes only, sampled at 'recordingSamplingRate' [nrTimePoints 1]
            % ixSpike = The index into spikeV that corresponds to simulated spikes.
            % neighborSpikeV = The spiking voltage observed in a
            % neighboring electtrode (same noise, same lfp, but a  different random draw of spikes)
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
                    nrSpk= o.duration*o.spikeRate;
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
                            ix =[ix;extra(take(1:nrToGenerate))]; %#ok<AGROW>
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
                plot(o.tRecord,o.vNeural/1e-6)
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
            % tacsShape    - sine,sawtooth, square
            % tdcsMean   -  Mean current [A]
            % tdcsRamp   - Duration of a ramp [s]
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

            % Determine the output
            switch upper(o.tacsShape)
                case 'SINE'
                    virtualOutputCurrent = o.tacsAmplitude*sin(2*pi*o.tacsFrequency*o.tTcs+o.tacsPhaseOffset);
                case 'SAWTOOTH'
                    % Sawtooth - with sine-phase (i.e. 0 at 0).
                    virtualOutputCurrent = o.tacsAmplitude*sawtooth(2*pi*o.tacsFrequency*o.tTcs+o.tacsPhaseOffset+pi);
                case 'SQUARE'
                    virtualOutputCurrent = o.tacsAmplitude*(sin(2*pi*o.tacsFrequency*o.tTcs+o.tacsPhaseOffset)>0);
                case 'NONE'
                    virtualOutputCurrent = zeros(1,numel(o.tTcs));
                otherwise
                    error('Unknown tACS shape %s',o.tacsShape);
            end

            if o.tdcsMean~=0
                virtualOutputCurrent = virtualOutputCurrent + o.tdcsMean*ones(numel(o.tTcs),1);
            end
            if o.tcsRamp >0
                inRamp = o.tTcs <= o.tcsRamp;
                nrRampSamples= find(inRamp,1,'last');
                scale = [linspace(0,1,nrRampSamples)'; ones(numel(o.tTcs)-2*nrRampSamples,1); linspace(1,0,nrRampSamples)'];
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

            if isa(o.zFun,'function_handle')
                % User-supplied function, pass the object
                z = o.zFun(o);
            elseif o.zAmplitude >0
                % Simulate a rhtyhm
                if o.tacsFrequency>0 && o.tacsFrequency/o.zFrequency == round(o.tacsFrequency/o.zFrequency)
                    fprintf(2,'Heartbeat is effectively locked to tACS - this is not a good simulation of a heartbeat... (use a frequency that is not a multiple of the tACS frequency)\n');
                end
                % Simulate a periodic "breathing/heartbeat" that changes the resistance - Pure sinusoid
                % Or with some long-term variability on a scale of 1 hour.
                frequency = o.zFrequency*(1+o.zVariability*sin(2*pi*(1/3600)*o.tSimulate));
                z = sin(2*pi*frequency.*o.tSimulate);
                if o.zDuration >0
                    % Create bumps with a given FWHM by raising the sine to
                    % a power of n
                    n = (-log(2)/log(cos(pi*o.zFrequency*o.zDuration)));
                    z(z<0) = 0;
                    z= z.^n;
                    z = z./max(z);
                end
                z = 1+ o.zAmplitude.*z;
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
            %simRecording - Simulate recording electrophysiogical signals produced by spikes, lfp, and tACS
            % Simulate the output of the amplifier that records the combination of a stimulation artifact
            % and an extracellular action potential trace.
            %
            % OPTIONAL EXTRA PARM/VALUE PAIRS
            % graph - Show a graph  [true]
            % showTime - limit graph to the first seconds [2]
            %
            % COMPUTES
            %  vRecord - The Voltage as recorded by the device. The first column represents
            %       the voltage during sham (including  additive noise ). The second column
            %       represents the same neural signal and additive noise,
            %       but now with the stimulation artifacts (including multiplicative noise).
            % vRecordNeighbor - Same as vRecord, but for a neihboring
            %                       electrode (Same lfp , different
            %                       spikes).
            %
            %  vRecordTacs- the applied stimulation voltage, as recorded by the
            %               recording device (and thus reflecting
            %               digitization of the stimulator).
            arguments
                o (1,1) artSim
                pv.graph (1,1) logical =false
                pv.showTime (1,1) double = 2
            end

            % Store the recorded vTcs assuming a linear amplifier, with a
            % highpass filter
            if o.highPass >0
                highPassFiltered  = highpass(o.vTcsActual,o.highPass,o.quasiContinuousSamplingRate);
            else
                highPassFiltered = o.vTcsActual;
            end
            o.vRecordTcs = decimate(highPassFiltered,o.quasiContinuousSamplingRate/o.recordingSamplingRate);

            %Add artifact, neural signal, and noise
            if o.pinkNoise
                additiveNoise= o.additive*coloredNoise(nrTimePoints=o.nrTimePoints,alpha=1);
            else
                additiveNoise = o.additive*randn(o.nrTimePoints,1);
            end
            multiplicativeNoise = 1+o.multiplicative.*randn(o.nrTimePoints,1);
            recV = nan(o.nrSamples,2);
            for i=1:2  % Electrode(2)  and neighbor (1)
                if i==1
                    simV = o.vNeuralNeighbor;
                else
                    simV = o.vNeural;
                end
                % Create two signals; 1 -> without tacs (.vTruth), 2-> with tACS
                simV = simV +additiveNoise + [zeros(size(simV)) o.vTcsActual.*multiplicativeNoise];

                % The recording hardware will filter this quasi continouous
                % signal before amplification/sampling. Here we model a
                % high pass hardware filter and an anti-aliasing filter (inside decimate).
                for j=1:size(simV,2)
                    if o.highPass >0
                        highPassFiltered  = highpass(simV(:,j),o.highPass,o.quasiContinuousSamplingRate);
                    else
                        highPassFiltered  = simV(:,j);
                    end
                    recV(:,j) = decimate(highPassFiltered,o.quasiContinuousSamplingRate/o.recordingSamplingRate);
                end
                % Amplification
                [amplifiedV] = artSim.amplify(recV,o.adcLinearRange,o.adcNonlinearity,o.adcGain);
                % Sampling  (assuming round to nearest bit)
                recV = o.adcVoltsPerBit*round(amplifiedV./o.adcVoltsPerBit);
                if i==1
                    o.vRecordNeighbor = recV;
                else
                    o.vRecord = recV;
                end
            end

            %% Summarize the results
            if pv.graph
                range = [prctile(amplifiedV,0.05) ; prctile(amplifiedV,99.5)]/1e-3;
                headroom = (2.^o.adcBitDepth-1) - 2*max(abs(bit));

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
                        [y,satPoint]= artSim.amplify(x,o.adcLinearRange,o.adcNonlinearity,o.adcGain,satLevel=99);
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
        function results = simPlotEEG(o,prms,arMode,pv)
            % Helper function to run the simulation, recording and artifact correction
            % and show the results for an EEG (or LFP) experiment. See
            % eegArtifacts.mlx for examples.
            arguments
                o (1,1) artSim
                prms (1,1) struct               % Artifact removal parameters
                arMode (1,:) string             = "FASTR"
                pv.freqRes (1,1) double         = 1
                pv.highCutOff (1,1) double      = 35 % Hz
                pv.lowCutOff (1,1) double       = 1 % Hz
                pv.tag (1,1) string             = ""
                pv.exportFig (1,1) logical      = false
                pv.errorThreshold (1,1) double  = 1 % More than 1% error is "forbidden".
                pv.fillHz (1,1) double          = 1;
                pv.axs  = []
                pv.tlim = [0 0.1] + mean(o.tRecord);
                pv.lineWidth = 2;
            end

            if isfield(prms,'tacsFrequency') && prms.tacsFrequency ~=o.tacsFrequency
                fprintf("Frequency in the simulation (%.2f) does not match the frequency in the ar (%.2f)\n",o.tacsFrequency,prms.tacsFrequency)
            end
            %% Initialize the simulation , simulate the recording, then perform artifact removal.
            o.rng               % Reset rng if .reproducible
            o.reset;            % Start fresh
            o.simLfp;           % 1. Generate LFP/EEG .
            o.simTCS;          % 2. Generate the TCS voltage
            o.simRecording;      % 3. Simulate recording
            % Now run artifact removal
            [vRecovered,results] = artifactRemoval(o.vContaminated,...
                prms, ...
                'groundTruth',o.vTruth,...
                'vTacsRecord',o.vRecordTcs,... % Include the tacs voltage in the PCA
                'recordingSamplingRate',o.recordingSamplingRate,...
                'mode',cellstr(arMode));

            %% Analyze the power spectrum
            [pwr,frequency] = pspectrum([o.vTruth vRecovered o.vContaminated ],o.recordingSamplingRate,'FrequencyLimits',[pv.lowCutOff pv.highCutOff],'FrequencyResolution',pv.freqRes);
            amplitude = sqrt(2)*sqrt(pwr); % Convert to muV amplitude.  (pspectrum is single sided PSD in rms,hence *sqrt(2))

            % Determine relative error in frequency domain
            hasSignal = amplitude(:,2)>1e-6;

            absError = amplitude(:,2)-amplitude(:,1);
            pctError = 100*(absError)./amplitude(:,1); % Divide by neural truth
            pctError(~hasSignal) = NaN;
            meanPctError = mean(abs(pctError),"omitmissing");
            results.amplitude = amplitude;
            results.pctError = pctError;
            results.frequency = frequency;
            results.meanPctError  = meanPctError;


            %% Assess the forbidden zones in the spectrum
            nrToFill = pv.fillHz/(results.frequency(3)-results.frequency(2)) ;
            [~,~,from, to] = fillRegions(abs(results.pctError)>pv.errorThreshold,nrToFill);
            forbidden = (results.frequency(to)-results.frequency(from));
            results.forbidden = sum(forbidden);
            if o.tacsFrequency==0
                inTacsBand = false;
            else
                if ~isempty(from) & ~isempty(to)
                    inTacsBand = from<find(results.frequency<o.tacsFrequency,1,"last") & to > find(results.frequency > o.tacsFrequency,1,"first");
                else
                    inTacsBand =false;
                end
            end
            if any(inTacsBand)
                results.forbiddenBand = forbidden(inTacsBand);
            else
                results.forbiddenBand = 0;
            end


            %% Generate a figure
            if isempty(pv.axs)
                [f,axs] = fig('Name',strjoin(arMode,'/')+pv.tag,'paperCols',1,'height',9,'nrRows',3,'byColumn',false);
            else
                axs = pv.axs;
            end
            if numel(axs)>=1 && isgraphics(axs(1))
                axes(axs(1)) % Time course
                yyaxis left  % Truth and recovered signal (small)
                tStay = o.tRecord >= pv.tlim(1)  &  o.tRecord < pv.tlim(2) ;
                h = plot(o.tRecord(tStay)-pv.tlim(1),[o.vTruth(tStay) vRecovered(tStay)]./1e-6,'LineWidth',1,'LineStyle','-');
                ylabel 'Amplitude (\muV)'
                ax =gca;
                ax.YColor ='k';
                hold on
                plot(xlim,zeros(1,2),'k-')
                ylim([-1 1]*max(abs(ylim)))
                yyaxis right % Signal with artifact

                h(3) = plot(o.tRecord(tStay)-pv.tlim(1),o.vContaminated(tStay)/1e-6,'LineWidth',1,'LineStyle','-');
                h(1).Color = 'r';
                h(2).Color = 'g';
                h(3).Color = 'b';
                h(1).LineWidth =pv.lineWidth;
                ylim([-1 1]*max(abs(ylim)))
                

                ax =gca;
                ax.YColor ='b';
                xlabel 'Time (s)'
                ylabel 'Amplitude (\muV)'

            end
            if  numel(axs)>=2 && isgraphics(axs(2))
                axes(axs(2)) % Spectrum
                yyaxis left  % Truth and recovered signal (small
                h = plot(frequency,amplitude(:,[1 2])./1e-6,'LineWidth',1,'LineStyle','-');
                ylabel 'Amplitude (\muV)'
                hold on
                plot(xlim,zeros(1,2),'k-')
                ax =gca;
                ax.YScale = 'Log';
                ax.YColor ='k';

                yyaxis right % Signal with artifact
                h(3) = plot(frequency,amplitude(:,3)/1e-6,'LineWidth',1,'LineStyle','-');
                h(1).Color = 'r';
                h(2).Color = 'g';
                h(3).Color = 'b';
                h(1).LineWidth =2;

                legend('Neural','Recovered','Recorded')

                ax =gca;
                ax.YScale = 'Log';
                ax.YColor ='b';
                set(gca,'XScale','Linear','XLim',[pv.lowCutOff pv.highCutOff]);
                xlabel 'Frequency (Hz)'
                ylabel 'Amplitude (\muV)'
            end
            if  numel(axs)>=3 && isgraphics(axs(3))
                axes(axs(3)); % Error as a percentage
                yyaxis left
                plot(frequency,pctError,'LineWidth',1,'Color','k');                
                xlabel 'Frequency (Hz)'
                ylabel 'Error (%)'
                set(gca,'XScale','Linear','YScale','Linear','XLim',[pv.lowCutOff pv.highCutOff],'YColor','k');
                ylim([-1 1]*max(abs(ylim)))
                hold on
                plot(xlim,zeros(1,2),'k-')
                yyaxis right
                plot(frequency,absError/1e-6,'LineWidth',1,'Color','m')
                set(gca,'YScale','Linear','YColor','m');                
                ylabel 'Error (\muV)'
               ylim([-1 1]*max(abs(ylim)))
            end

            fprintf('%s (r=%.2f rmse = %.0f/%.0f muV ae = %.1f%%)\n',strjoin(arMode,'/'),results.r.(arMode),results.mse.(arMode),o.additive/1e-6,meanPctError);

            % Export figure for paper
            if pv.exportFig
                name = sprintf('figure%s.png',pv.tag);
                exportgraphics(f,fullfile('../docs/',name),"Colorspace",'rgb','ContentType','vector','Resolution',300);
            end
        end


        function [sortResults,falsePositives,falseNegatives,spikes,vRecovered] = simPlotSpikes(o,prms,arMode,pv)
            % Helper function to run the simulation, recording and artifact correction
            % and show the results for an single unit recording experiment. See
            % spikeArtifacts.mlx for examples.
            arguments
                o (1,1) artSim
                prms (1,1) struct       % Artifact removal parms
                arMode (1,1) string = "FASTR"
                pv.lowCut (1,1) double = 300
                pv.highCut (1,1) double = 3000 % Hz
                pv.filterOrder (1,1) double = 1
                pv.threshold (1,:) double =1
                pv.sharedThreshold (1,1) logical= true
                pv.spikeSortingMode (1,1) string = "MATCHSHAM"
                pv.nrPhaseBins (1,1) double = 8
                pv.phaseFromRaw (1,1) logical = true;
                pv.showMissed (1,1) logical  = false
                pv.tag (1,1) string = ""
                pv.exportFig (1,1) logical = false
                pv.axs  = []
                pv.xlim  = 4+[0 0.250];
            end
            colors = 'br';
            o.rng               % Reset rng if .reproducible
            o.reset;            % Start clean
            o.simLfp;           % 1. Generate LFP/EEG .
            o.simSpikes;        % 2. Generate spikes
            o.simTCS;          % 3. Generate the tACS voltage
            o.simRecording;      % 4. Simulate recording


            %% Run artifact removal
            [vRecovered] = artifactRemoval(o.vContaminated,...
                prms, ...
                'tacsFrequency',o.tacsFrequency,...
                'groundTruth',o.vTruth,...
                'vTacsRecord',o.vRecordTcs,... % Include the tacs voltage in the PCA
                'recordingSamplingRate',o.recordingSamplingRate,...
                'mode',{arMode});

            %% Sort spikes
            [sortResults,spikes] =sortSpikes(o,vRecovered,...
                'sortMode',pv.spikeSortingMode, ...
                'spikeBand',[pv.lowCut pv.highCut], ...
                'filterOrder',pv.filterOrder, ...
                'threshold',pv.threshold, ...
                'sharedThreshold',pv.sharedThreshold, ...
                'phaseFromRaw',pv.phaseFromRaw, ...
                'nrPhaseBins',pv.nrPhaseBins);


            %% Visualize
            % Generate a figure
            if isempty(pv.axs)
                [f,axs] = fig('Name',"Spikes" + arMode,'paperCols',2,'height',9,'nrRows',2,'nrCols',2,'byColumn',true);
            else
                axs = pv.axs;
            end

            % A/B - Filtered voltage and spike detection
            for i=1:2
                if isgraphics(axs(i))
                    axes(axs(i)); %#ok<LAXES>
                    h0= plot(o.tRecord,sortResults(i).voltage/(o.adcGain*1e-6),'LineWidth',0.5);
                    h1 = plot(o.tSimulate(o.ixSpike),zeros(size(o.ixSpike)),'go','MarkerSize',7); % Ground truth
                    h2= plot(o.tRecord(sortResults(i).ixDetected),zeros(size(sortResults(i).ixDetected)),'g*','MarkerSize',7);
                    ylim(o.spikePeak/1e-6*[-1 1])
                    plot(pv.xlim,spikes.info.detect.thresh/(o.adcGain*1e-6)*[1 1],'k')
                    xlabel 'Time (s)'
                    ylabel ('Filtered Voltage (\muV)')
                    if i==1
                        legend([h0 h1 h2],'V','Spk','Detect','Orientation','horizontal','Location','north')
                    else
                        h3= plot(o.tTcs,max(ylim)*o.vTcs/max(o.vTcs),'b:','LineWidth',1);
                        legend (h3,'tACS')
                    end
                    xlim(axs(i),pv.xlim)
                end
            end
            % C  - Polar spike histogram
            if numel(axs)>2 && isgraphics(axs(3))
                axs(3).Visible = 'off';
                axs(3) = polaraxes('Position',axs(3).Position);

                h = nan(1,2);
                maxRho =0;
                for i=1:2
                    theta   = [sortResults(i).binCenters;sortResults(i).binCenters(1)];
                    rho     = [sortResults(i).spikeCount;sortResults(i).spikeCount(1)]./[sortResults(i).phaseCount;sortResults(i).phaseCount(1)];
                    h(i)= polarplot(axs(3),theta,rho,'Color',colors(i),'LineStyle','-','LineWidth',1);
                    hold on
                    polarplot(axs(3),theta,rho, '.','Color',colors(i));
                    maxRho = max([maxRho;rho]);
                end
                axs(3).ThetaTick = [0 45 75 90 135 180 225 270 315];
                axs(3).ThetaTickLabel = {'0','45','','90','135','180','225','270','315'};
                axs(3).RLim = [0 maxRho];
                axs(3).RTick = [0 maxRho];
                axs(3).RTickLabel = [];
                axs(3).RAxisLocation =75;
                %axs(3).RAxis.Label.String ='p(spike)';
                %axs(3).RAxis.Label.Rotation = 70;
                %axs(3).RAxis.Label.Position= [70 0.5];

                %axs(3).ThetaAxis.Label.String ='tACS phase';
                %axs(3).ThetaAxis.Label.Position =[0 1.1];
                %axs(3).ThetaAxis.Label.Rotation   = 90;
                legend(h, 'Sham','tACS','Location','SouthEastOutside')
            end

            % D - waveforms for each phase
            if numel(axs)>3 && isgraphics(axs(4))
                axes(axs(4))
                axs(4).Visible = 'off';
                if pv.showMissed
                    m =   cat(3,sortResults(1).meanMissedWF(25:45,:),sortResults(2).meanMissedWF(25:45,:));
                    sd = cat(3,sortResults(1).sdMissedWF(25:45,:),sortResults(2).sdMissedWF(25:45,:));
                else
                    m =   cat(3,sortResults(1).meanWF(25:45,:),sortResults(2).meanWF(25:45,:));
                    sd = cat(3,sortResults(1).sdWF(25:45,:),sortResults(2).sdWF(25:45,:));
                end

                plotWaveformPhase(sortResults(1).binCenters, m,sd,...
                    'ax',axs(4),'UseBins',1:2:pv.nrPhaseBins, ...
                    'radius',0.25,'width',0.4,'height',0.4, ...
                    'colors',colors);
            end

            if pv.exportFig
                name = sprintf('figure%s.png',pv.tag);
                exportgraphics(f,fullfile('../docs/',name),"Colorspace",'rgb','ContentType','vector','Resolution',300);
            end

            % Summarize results on the command line
            trueNrSpikes = numel(sortResults(1).found); % Same for both
            falsePositives = 100*[sum(sortResults(1).invented) ;sum(sortResults(2).invented)]/trueNrSpikes;
            falseNegatives = 100*[sum(~sortResults(1).found) ;sum(~sortResults(2).found)]/trueNrSpikes;
            fprintf('False Positive. Sham: %.0f%%  tACS: %.0f%% \n',falsePositives(1),falsePositives(2))
            fprintf('False Negative. Sham: %.0f%%  tACS: %.0f%% \n',falseNegatives(1),falseNegatives(2))
            fprintf('PLV. Sham: %.2f  tACS: %.2f \n',sqrt(max(0,sortResults(1).ppc)),sqrt(max(0,sortResults(2).ppc)))
            fprintf('PLV2. Sham: %.2f  tACS: %.2f \n',sqrt(max(0,sortResults(1).ppc2)),sqrt(max(0,sortResults(2).ppc2)))

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
                pv.sharedThreshold (1,1) logical = false
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
                spikes = ss_detect({filteredVoltage(:,1)},spikes);

                %Trial 2 = recovered signal
                spikes.params.thresh =  std(filteredVoltage(:,2))* -pv.threshold;
                spikes.params.detect_method ='manual';
                spikes = ss_detect({filteredVoltage(:,2)},spikes);
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
                results(i).ixDetected = ixDetectedSpike(isAssignedSpike);
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
                voltageInBand = downsample(o.vRecordNeighbor,o.recordingSamplingRate/1000);
            else
                % Estimate phase from the tACS signal.
                voltageInBand = downsample(o.vRecordTcs,o.recordingSamplingRate/1000);
                voltageInBand = repmat(voltageInBand,[1 2]);
            end
            [B,A] = butter(3,(o.tacsFrequency+[-1 +1])/(1000/2),'bandpass');
            voltageInBand = filtfilt(B,A,voltageInBand );
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