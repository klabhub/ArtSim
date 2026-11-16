classdef defaultArtifactRemovalTest < matlab.unittest.TestCase
    % Characterization tests for the default EEG-demo artifact removers.
    %
    % The one-minute baseline was generated from the settings in demo/eeg.m,
    % using the current artSim property name amplifierLinearRange.

    properties
        baseline
        actualFbar
        actualRbar
        actualFastr
    end

    methods (TestClassSetup)
        function createDefaultDemoResults(testCase)
            rootFolder = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                rootFolder,IncludingSubfolders=true));
            testCase.baseline = load(fullfile(rootFolder,'test','baselines', ...
                'defaultArtifactRemovalBaseline.mat'));

            defaults = [];
            defaults.reproducible = true;
            defaults.stimulatorSamplingRate = 1e3;
            defaults.currentResolution = 1e-6;
            defaults.recordingSamplingRate = 5e3;
            defaults.amplifierLinearRange = 200e-3;
            defaults.adcNonlinearity = 100;
            defaults.adcRange = 200e-3;
            defaults.adcBitDepth = 24;
            defaults.amplifierGain = 1;
            defaults.highPass = 0.01;
            defaults.lowPass = 100;
            defaults.duration = 60;
            defaults.tacsAmplitude = 1e-3;
            defaults.resistance = 15;
            defaults.tacsFrequency = 10;
            defaults.additive = 100e-6;
            defaults.lfpFrequency = 24;
            defaults.lfpAmplitude = 10e-6;
            defaults.lfpPeakWidth = 1;

            o = artSim(defaults);
            o.sim;

            fastr = [];
            fastr.name = "FASTR";
            fastr.options.slack = 5e-3;
            fastr.options.nrSegsPerWindow = 10;
            fastr.options.pcaNr = 2;
            fastr.options.pcaVarExplained = Inf;
            fastr.options.pcaNrSegsPerWindow = Inf;
            fastr.options.referenceSegment = 1;

            fbar = [];
            fbar.name = "FBAR";
            fbar.options.frequencyLimits = [-1 1];

            rbar = [];
            rbar.name = "RBAR";
            rbar.options = [];

            testCase.actualFbar = artifactRemoval( ...
                o.vContaminated,fbar, ...
                recordingSamplingRate=defaults.recordingSamplingRate, ...
                segmentDuration=3,tacsFrequency=defaults.tacsFrequency).signal;
            testCase.actualRbar = artifactRemoval( ...
                o.vContaminated,rbar, ...
                recordingSamplingRate=defaults.recordingSamplingRate, ...
                segmentDuration=3,tacsFrequency=defaults.tacsFrequency, ...
                vRecordTcs=o.vRecordTcs).signal;
            testCase.actualFastr = artifactRemoval( ...
                o.vContaminated,fastr, ...
                recordingSamplingRate=defaults.recordingSamplingRate, ...
                segmentDuration=3,tacsFrequency=defaults.tacsFrequency).signal;
        end
    end

    methods (Test)
        function testFbarMatchesCurrentBaseline(testCase)
            testCase.verifyEqual(testCase.actualFbar,testCase.baseline.fbar, ...
                AbsTol=1e-7);
        end

        function testRbarMatchesCurrentBaseline(testCase)
            testCase.verifyEqual(testCase.actualRbar,testCase.baseline.rbar, ...
                AbsTol=1e-7);
        end

        function testFastrMatchesCurrentBaseline(testCase)
            testCase.verifyEqual(testCase.actualFastr,testCase.baseline.fastr, ...
                AbsTol=1e-7);
        end
    end
end
