% runAllTests  Run the complete ArtSim test suite.

testFolder = fileparts(mfilename('fullpath'));
results = runtests(testFolder);

passed = [results.Passed];
failed = [results.Failed];
incomplete = [results.Incomplete];

fprintf('Tests: %d total, %d passed, %d failed, %d incomplete.\n', ...
    numel(results),sum(passed),sum(failed),sum(incomplete));

if any(failed | incomplete)
    error('ArtSim:TestsFailed', ...
        'The test suite contains failed or incomplete tests.');
end
