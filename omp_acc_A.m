clear;

% File parameters
mName = mfilename('fullpath');
fileName = get_name(mName);

% Default parameters 
lambda = 0.5;                   % wavelength (cm)
size_of_hand = 16;              % hand box size (cm)
scatterer_grid_size = 32;       % number of scatterers in each axis 
grid_height = 30;               % grid height (cm)
num_antennas = 46;              % number of antennas in array (on each side)
random_seed = 'shuffle';        % 'shuffle' for random; any int for reproducibility
nIter = 1000;                   % number of attempts 

% Random generator
rStr = RandStream('mcg16807', 'Seed', random_seed);

% Generate A
[A, numTX, numRX] = gen_A(lambda, size_of_hand, scatterer_grid_size, grid_height, num_antennas);
iF = gen_iF(scatterer_grid_size);
AiF = A * iF;

N = size(AiF, 2);
s = 102;
step = 2;
start = 8;
stop = 32;

accuracy = zeros(N, 1);
for Neff = start:step:stop

    right = 0;
    M = Neff ^ 2;

    for i = 1:nIter
        % Progress
        fprintf('Iter: %d, M = %d\n', i, M);

        txPerm = randperm(rStr, 2 * num_antennas, Neff)';
        rxPerm = randperm(rStr, 2 * num_antennas, Neff)';
        ATrim = trim_A(AiF, txPerm, rxPerm, numRX);

        % Generate sparse x
        [x, ~] = sparse_x(rStr, N, s, true);
    
        % Sample 
        y = ATrim * x;
    
        % Perform OMP
        [xHat, ~] = omp(y, ATrim);
        error = sum(abs(xHat - x) .^ 2);
    
        if error < 1e-10
            right = right + 1;
        end
    end

    accuracy(M) = right / nIter;
    save(fileName, 'nIter', 'N', 's', 'accuracy', 'step', 'start', 'stop');

end
