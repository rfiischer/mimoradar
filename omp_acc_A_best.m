clear;

% File parameters
mName = mfilename('fullpath');
fileName = get_name(mName);

% Default parameters
random_seed = 'shuffle';        % 'shuffle' for random; any int for reproducibility
nIter = 10;                   % number of attempts 
post_process = @fill_gap;       % post processing function
experimentIdx = 6;              % experiment identifier 
description = 'default';        % simulation description 

% Simulation parameters 
s = 0.1;
step = 2;
start = 24;
stop = 26;

% Random generator
rStr = RandStream('mcg16807', 'Seed', random_seed);

% 'A' database directory
baseName = sprintf('A_database_%d//', experimentIdx);
if ~exist(baseName, 'dir')
    error('No database directory with this number %d', experimentIdx);
end

valueRange = start:step:stop;
accuracy = zeros(size(valueRange, 2), 2);
for m = 1:length(valueRange)

    right1 = 0;
    right2 = 0;
    Neff = valueRange(m);
    M = Neff ^ 2;
    
    fileNameA = fullfile(baseName, sprintf('A_%d.mat', Neff));
    if ~exist(fileNameA, 'file')
        error('Matrix not found: A_%d', Neff);

    else
        loadA = load(fileNameA);
        ACand = loadA.ACand;
        AMean = loadA.AMean;
        if isinteger(ACand) && isinteger(AMean)
            error('Either ACand or AMean is zero.');
        else
            N = size(ACand, 2);
        end

    end
    
    for i = 1:nIter
        % Progress
        fprintf('Iter: %d, M = %d\n', i, M);

        % Generate sparse x
        [x, ~] = sparse_x(rStr, N, s, true);
    
        % Sample 
        y1 = ACand * x;
        y2 = AMean * x;
    
        % Perform OMP
        [xHat1, ~] = omp(y1, ACand, post_process);
        [xHat2, ~] = omp(y2, AMean, post_process);
        error1 = sum(abs(xHat1 - x) .^ 2);
        error2 = sum(abs(xHat2 - x) .^ 2);
    
        if error1 < 1e-10
            right1 = right1 + 1;
        end

        if error2 < 1e-10
            right2 = right2 + 1;
        end
    end

    accuracy(m, 1) = right1 / nIter;
    accuracy(m, 2) = right2 / nIter;
    save(fileName, 'nIter', 'N', 's', 'accuracy', 'step', 'start', 'stop', ...
        'description', 'experimentIdx', 'valueRange');

    

end
