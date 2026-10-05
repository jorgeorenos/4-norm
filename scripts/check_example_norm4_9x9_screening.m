% Exercise a small SO(9) histogram sweep with the screened settings.
scriptDirectory = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDirectory);
oldVisibility = get(groot, 'DefaultFigureVisible');
set(groot, 'DefaultFigureVisible', 'off');
try
    source = fileread(fullfile(scriptDirectory, 'example_norm4_9x9.m'));
    source = strrep(source, "scriptDirectory = fileparts(mfilename('fullpath'));", ...
        "scriptDirectory = fileparts(which('check_example_norm4_9x9_screening'));");
    source = strrep(source, 'numMatrices = 1000;', 'numMatrices = 50;');
    source = strrep(source, "'MaxIterations', 1000", "'MaxIterations', 30");
    eval(source);
    assert(isequal(size(norm4_norms.starts1), [50 1]));
    assert(isequal(size(norm4_norms.starts10), [50 1]));
    assert(isequal(size(norm4_norms.starts25), [50 1]));
    assert(isequal(size(norm4_norms.starts50), [50 1]));
    assert(isequal(size(norm4_norms.starts100), [50 1]));
    assert(options.ScreenIterations == 4 && options.NumFinalists == 3);
    assert(numel(findobj(gcf, 'Type', 'histogram')) == 9);
    rng(42, 'twister');
    referenceQ = zeros(9, 9, 50);
    for page = 1:50
        referenceQ(:,:,page) = haar_so(9);
    end
    assert(isequal(referenceQ, rotations));
    for count = startCounts
        expected = struct('NumRandomStarts', count, 'MaxIterations', 30, ...
            'MaxWorkingMemoryMB', 128);
        if count >= 50
            expected.ScreenIterations = 4;
            expected.NumFinalists = 3;
        end
        reference = compute_induced_norm(referenceQ, @norm_4, expected);
        assert(max(abs(reference - norm4_norms.(sprintf('starts%d', count)))) < 1e-12);
    end
    close all;
catch err
    close all;
    set(groot, 'DefaultFigureVisible', oldVisibility);
    rethrow(err);
end
set(groot, 'DefaultFigureVisible', oldVisibility);
disp('SO(9) screened histograms: OK');
