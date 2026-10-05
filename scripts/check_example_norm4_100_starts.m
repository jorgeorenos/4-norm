% Verify the focused SO(9) experiment and its reported workload.
scriptDirectory = fileparts(mfilename('fullpath'));
oldFigures = numel(findall(groot, 'Type', 'figure'));
report = evalc('run(fullfile(scriptDirectory, ''example_norm4_100_starts.m''));');
assert(isequal(size(estimates), [1000, 1]));
assert(all(isfinite(estimates) & estimates >= 1-1e-12 & ...
    estimates <= 9^(1/4)+1e-12));
assert(numMatrices == 1000 && numStarts == 100 && elapsedSeconds > 0);
assert(contains(report, '1000') && contains(report, '100') && ...
    contains(report, 'seconds'));
assert(numel(findall(groot, 'Type', 'figure')) == oldFigures);
disp('Focused SO(9) example: OK');
