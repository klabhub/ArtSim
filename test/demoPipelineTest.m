classdef demoPipelineTest < matlab.unittest.TestCase
    % Smoke-tests the reusable EEG analysis function called by eeg.m.

    methods (TestClassSetup)
        function addDemoFolderToPath(testCase)
            rootFolder = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                rootFolder,IncludingSubfolders=true));
        end
    end

    methods (Test)
        function testAnalyzeEegAcceptsArtifactRemovalResult(testCase)
            o = artSim;
            o.duration = 10;
            o.recordingSamplingRate = 1000;
            o.stimulatorSamplingRate = 100;
            o.reproducible = true;
            o.tacsFrequency = 10;
            o.tacsAmplitude = 1e-3;
            o.resistance = 1;
            o.sim;

            artifactResult = artifactRemoval(o.vContaminated, ...
                struct('name','NOP'), ...
                recordingSamplingRate=o.recordingSamplingRate, ...
                segmentDuration=o.duration,groundTruth=o.vTruth);
            analysis = analyzeEEG(o,artifactResult, ...
                freqRes=1,lowCutOff=1,highCutOff=35,axs=[],exportFig=false);

            testCase.verifyTrue(isfield(analysis,'frequency'));
            testCase.verifyTrue(isfield(analysis,'pctError'));
            testCase.verifyEqual(numel(analysis.frequency), ...
                size(analysis.amplitude,1));
            testCase.verifyTrue(all(isfinite(analysis.frequency)));
        end
    end
end
