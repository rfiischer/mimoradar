clear;

% Default parameters 
lambda = 0.5;                   % wavelength (cm)
scatterer_grid_size = 32;       % number of scatterers in each axis 
scatterer_spacing = 0.5;        % scatterer spacing (cm)
grid_height = 30;               % grid height (cm)
num_antennas = 46;              % number of antennas in array (on each side)
antenna_spacing = lambda / 2;   % antenna spacing (cm)
num_eff_antennas = 20;          % effective number of antennas (both sides) chosen randomly 
random_seed = 'shuffle';        % 'shuffle' for random; any int for reproducibility
numIter = 1000;                 % number of attempts 

% Random generator
rStr = RandStream('mcg16807', 'Seed', random_seed);

% File parameters
mName = mfilename('fullpath');
fileName = get_name(mName);

% Build A and iF
[txPoints, rxPoints, scattererPoints] = generate_setup(antenna_spacing, ...
    num_antennas, scatterer_spacing, scatterer_grid_size, grid_height);
numRX = size(rxPoints, 1);
A = gen_A(lambda, txPoints, rxPoints, scattererPoints);
iF = gen_iF(scatterer_grid_size);
AiF = A * iF;

[ACand, AMean, co, minCand, meanCand, txPerm, rxPerm, i] = random_search_A(rStr, AiF, numRX, ...
    num_antennas, num_eff_antennas, numIter, fileName, true);

save(fileName, 'lambda', 'size_of_hand', 'scatterer_grid_size', 'grid_height', ...
            'num_antennas', 'num_eff_antennas', 'random_seed', 'numIter', 'co', 'minCand', 'ACand', ...
            'AMean', 'meanCand', 'txPerm', 'rxPerm', 'scattererPoints', 'txPoints', 'rxPoints', 'i');
