clear;

% File parameters
mName = mfilename('fullpath');
fileName = get_name(mName);

% Default parameters 
lambda = 0.5;                   % wavelength (cm)
scatterer_grid_size = 32;       % number of scatterers in each axis 
scatterer_spacing = 0.5;        % scatterer spacing (cm)
grid_height = 30;               % grid height (cm)
num_antennas = 46;              % number of antennas in array (on each side)
antenna_spacing = lambda / 2;   % antenna spacing (cm)
random_seed = 'shuffle';        % 'shuffle' for random; any int for reproducibility
nIter = 1000;                   % number of attempts 
post_process = @fill_gap;       % post processing function
description = 'default';        % simulation description 

% Simulation parameters
s = 0.1;
step = 2;
start = 8;
stop = 32;

% Random generator
rStr = RandStream('mcg16807', 'Seed', random_seed);

% Generate A
[txPoints, rxPoints, scattererPoints] = generate_setup(antenna_spacing, ...
    num_antennas, scatterer_spacing, scatterer_grid_size, grid_height);
numRX = size(rxPoints, 1);
A = gen_A(lambda, txPoints, rxPoints, scattererPoints);
iF = gen_iF(scatterer_grid_size);
AiF = A * iF;

N = size(AiF, 2);
valueRange = start:step:stop;
accuracy = zeros(size(valueRange, 2), 1);
for m = 1:length(valueRange)

    right1 = 0;
    right2 = 0;
    Neff = valueRange(m);
    M = Neff ^ 2;

    for i = 1:nIter
        % Progress
        fprintf('Iter: %d, M = %d\n', i, M);

        txPerm = randperm(rStr, 2 * num_antennas, Neff);
        rxPerm = randperm(rStr, 2 * num_antennas, Neff);
        ATrim = trim_A(AiF, txPerm, rxPerm, numRX);

        % Generate sparse x
        [x, ~] = sparse_x(rStr, N, s, true);
    
        % Sample 
        y = ATrim * x;
    
        % Perform OMP
        [xHat, ~] = omp(y, ATrim, post_process);
        error = sum(abs(xHat - x) .^ 2);
    
        if error < 1e-10
            right1 = right1 + 1;
        end
    end

    accuracy(m) = right1 / nIter;
    save(fileName, 'lambda', 'scatterer_grid_size', 'scatterer_spacing', ...
        'grid_height', 'num_antennas', 'antenna_spacing', 'random_seed', ...
        'nIter', 'post_process', 'description', 'start', 'step', 'stop', 'valueRange', ...
        'accuracy');


end
