classdef signalUtilitiesTest < matlab.unittest.TestCase
    % Tests signal utilities used by the EEG and spike analysis paths.

    methods (TestClassSetup)
        function addSourceFolder(testCase)
            rootFolder = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                rootFolder,IncludingSubfolders=true));
        end
    end

    methods (Test)
        function testFftRoundTrip(testCase)
            t = (0:999)'/1000;
            x = sin(2*pi*7*t) + 0.2*cos(2*pi*31*t);
            [ft,freq] = fftReal(x,1000);
            [recovered,time] = ifftReal(ft,freq,1000);

            testCase.verifySize(ft,[501 1]);
            testCase.verifySize(time,[1 1000]);
            testCase.verifyEqual(recovered,x,AbsTol=1e-12);
        end

        function testColoredNoiseHasRequestedShapeAndNormalization(testCase)
            rng(17);
            [signal,time] = coloredNoise( ...
                nrTimePoints=1000,samplingRate=1000,alpha=1);

            testCase.verifySize(signal,[1000 1]);
            testCase.verifySize(time,[1 1000]);
            testCase.verifyEqual(mean(signal),0,AbsTol=1e-12);
            testCase.verifyEqual(std(signal),1,AbsTol=1e-12);
        end

        function testFillRegionsFillsShortGaps(testCase)
            input = logical([1 0 0 1 0 1]);
            [filled,nrInRegion,startIx,stopIx] = fillRegions(input,2);

            testCase.verifyEqual(filled,true(6,1));
            testCase.verifyEqual(nrInRegion,6);
            testCase.verifyEqual(startIx,1);
            testCase.verifyEqual(stopIx,6);
        end

        function testQcIsPerfectForIdenticalSignals(testCase)
            x = sin(2*pi*(0:99)'/100);
            [correlation,mse] = qc(x,x);

            testCase.verifyEqual(correlation,1,AbsTol=1e-12);
            testCase.verifyEqual(mse,0,AbsTol=1e-12);
        end
    end
end
