clear;

% Default parameters 
lambda = 0.5;                   % wavelength (cm)
size_of_hand = 16;              % hand box size (cm)
scatterer_grid_size = 32;       % number of scatterers in each axis 
grid_height = 30;               % grid height (cm)
num_antennas = 46;              % number of antennas in array (on each side)
num_eff_antennas = 20;          % effective number of antennas (both sides) chosen randomly 
random_seed = 'shuffle';        % 'shuffle' for random; any int for reproducibility
numIter = 1000;                 % number of attempts 

% Random generator
rStr = RandStream('mcg16807', 'Seed', random_seed);

% File parameters
mName = mfilename('fullpath');
fileName = get_name(mName);

% Build A and iF
[A, numTX, numRX, scattererPoints, txPoints, rxPoints] = gen_A(lambda, size_of_hand, ...
    scatterer_grid_size, grid_height, num_antennas);
iF = gen_iF(scatterer_grid_size);
AiF = A * iF;
N = size(A, 2);

co = zeros(numIter, 1);
minCand = 1;
ACand = 0;
for i = 1:numIter
    fprintf('i: %d\n', i);

    % Pick random antennas
    txPerm = randperm(rStr, 2 * num_antennas, num_eff_antennas)';
    rxPerm = randperm(rStr, 2 * num_antennas, num_eff_antennas)';

    ATrim = trim_A(AiF, txPerm, rxPerm, numRX);

    co(i) = coherence(ATrim);

    if co(i) < minCand
        minCand = co(i);
        ACand = ATrim;
        save(fileName, "lambda", "size_of_hand", "scatterer_grid_size", "grid_height", ...
            "num_antennas", "num_eff_antennas", "random_seed", "numIter", "co", "minCand", "ACand", ...
            "txPerm", "rxPerm", "scattererPoints", "txPoints", "rxPoints", "i");
    end
    
end

save(fileName, "lambda", "size_of_hand", "scatterer_grid_size", "grid_height", ...
            "num_antennas", "num_eff_antennas", "random_seed", "numIter", "co", "minCand", "ACand", ...
            "txPerm", "rxPerm", "scattererPoints", "txPoints", "rxPoints", "i");
