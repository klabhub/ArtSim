classdef artifactRemovalTest < matlab.unittest.TestCase
    % Tests artifact-removal behavior used by both demo Live Scripts.

    methods (TestClassSetup)
        function addSourceFolder(testCase)
            rootFolder = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                rootFolder,IncludingSubfolders=true));
        end
    end

    methods (Test)
        function testNopPreservesSignalAndReportsSteps(testCase)
            x = (1:1000)'/1000;
            step = struct('name','NOP');
            result = artifactRemoval(x,step, ...
                recordingSamplingRate=1000,segmentDuration=1,groundTruth=x);

            testCase.verifyEqual(result.signal,x,AbsTol=0);
            testCase.verifyTrue(isfield(result,'step'));
            testCase.verifyTrue(isfield(result,'common'));
            testCase.verifySize(result.signal,size(x));
        end

        function testFbarPreservesLengthAndFiniteValues(testCase)
            t = (0:999)'/1000;
            x = 2e-3*sin(2*pi*10*t) + 1e-4*sin(2*pi*3*t);
            step = struct('name','FBAR', ...
                'options',struct('frequencyLimits',[-1 1]));
            result = artifactRemoval(x,step, ...
                recordingSamplingRate=1000,segmentDuration=1);

            testCase.verifySize(result.signal,size(x));
            testCase.verifyTrue(all(isfinite(result.signal)));
        end
    end
end
