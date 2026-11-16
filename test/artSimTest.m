classdef artSimTest < matlab.unittest.TestCase
    % Tests the simulation paths used by the EEG and spike demos.

    methods (TestClassSetup)
        function addSourceFolder(testCase)
            rootFolder = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                rootFolder,IncludingSubfolders=true));
        end
    end

    methods (Test)
        function testEegSimulationProducesAlignedSignals(testCase)
            o = artSim;
            o.duration = 1;
            o.recordingSamplingRate = 1000;
            o.stimulatorSamplingRate = 100;
            o.reproducible = true;
            o.tacsFrequency = 10;
            o.tacsAmplitude = 1e-3;
            o.resistance = 1;
            o.lfpAmplitude = 100e-6;
            o.lfpFrequency = 8;
            o.sim;

            testCase.verifySize(o.vTruth,[1000 1]);
            testCase.verifySize(o.vContaminated,[1000 1]);
            testCase.verifySize(o.tRecord,[1000 1]);
            testCase.verifyTrue(all(isfinite(o.vTruth)));
            testCase.verifyTrue(all(isfinite(o.vContaminated)));
        end

        function testSpikeSimulationProducesSpikes(testCase)
            o = artSim;
            o.duration = 2;
            o.recordingSamplingRate = 1000;
            o.stimulatorSamplingRate = 100;
            o.reproducible = true;
            o.spikePeak = 1e-3;
            o.spikeRate = 10;
            o.sim;

            testCase.verifyGreaterThan(o.nrSpikes,0);
            testCase.verifySize(o.vNeural,[200000 1]);
            testCase.verifySize(o.vTruth,[2000 1]);
            testCase.verifyEqual(numel(o.ixSpike),o.nrSpikes);
        end
    end
end
