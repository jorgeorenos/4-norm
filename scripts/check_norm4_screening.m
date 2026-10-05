% Screening must preserve the prescribed starts and approximate the full run.
scriptDirectory = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDirectory);
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
full = compute_induced_norm(Q, @norm_4, base);
endState = rng;
screened = base;
screened.ScreenIterations = 4;
screened.NumFinalists = 3;
rng(state);
fast = compute_induced_norm(Q, @norm_4, screened);
assert(isequal(size(fast),[20 1]) && max(abs(fast-full))<1e-10);
assert(isequal(rng,endState));
rng(state);
equivalent = compute_induced_norm(Q, @(x) norm_4(x), screened);
assert(max(abs(equivalent-fast)) < 1e-12 && isequal(rng,endState));
assert_error(@() compute_induced_norm(Q(:,:,1), @norm_4, ...
    struct('ScreenIterations',-1)), 'compute_induced_norm:InvalidCount');
assert_error(@() compute_induced_norm(Q(:,:,1), @norm_4, ...
    struct('NumRandomStarts',2,'ScreenIterations',4,'NumFinalists',3)), ...
    'compute_induced_norm:InvalidCount');
assert_error(@() compute_induced_norm(Q(:,:,1), @norm_4, ...
    struct('NumFinalists',0)), 'compute_induced_norm:InvalidCount');
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
