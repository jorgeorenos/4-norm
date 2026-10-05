% Screening must preserve the prescribed starts and approximate the full run.
scriptDirectory = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptDirectory));
addpath(fullfile(projectRoot, 'src', '4-norm'));
addpath(fullfile(projectRoot, 'src', 'helpers'));
rng(71, 'twister');
Q = zeros(9,9,20);
for page = 1:20
    Q(:,:,page) = haar_so(9);
end
state = rng;
base = struct('NumRandomStarts',100,'MaxIterations',1000, ...
    'MaxWorkingMemoryMB',0.6);
rng(state);
full = compute_4_norm(Q, base);
endState = rng;
screened = base;
screened.ScreenIterations = 4;
screened.NumFinalists = 3;
rng(state);
fast = compute_4_norm(Q, screened);
assert(isequal(size(fast),[20 1]) && max(abs(fast-full))<1e-10);
assert(isequal(rng,endState));
assert_error(@() compute_4_norm(Q(:,:,1), ...
    struct('ScreenIterations',-1)), 'compute_4_norm:InvalidCount');
assert_error(@() compute_4_norm(Q(:,:,1), ...
    struct('NumRandomStarts',2,'ScreenIterations',4,'NumFinalists',3)), ...
    'compute_4_norm:InvalidCount');
assert_error(@() compute_4_norm(Q(:,:,1), ...
    struct('NumFinalists',0)), 'compute_4_norm:InvalidCount');
disp('Screened induced 4-norm checks: OK');

function assert_error(call, identifier)
    try
        call();
    catch err
        assert(strcmp(err.identifier, identifier), err.message);
        return;
    end
    error('check_norm4_screening:ExpectedError', 'Expected %s.', identifier);
end
