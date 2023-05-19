clear;

% Parameters 
lambda = 0.5;                   % wavelength (cm)
scatterer_grid_size = 32;       % number of scatterers in each axis 
scatterer_spacing = 0.5;        % scatterer spacing (cm)
grid_height = 30;               % grid height (cm)
num_antennas = 46;              % number of antennas in array (on each side)
antenna_spacing = lambda / 2;   % antenna spacing (cm)
num_eff_antennas = 20;          % effective number of antennas (both sides) chosen randomly 
random_seed = 0;                % 'shuffle' for random; any int for reproducibility

% Random generator
rStr = RandStream('mcg16807', 'Seed', random_seed);

% Get geometry points
[txPoints, rxPoints, scattererPoints] = generate_setup(antenna_spacing, num_antennas, ...
    scatterer_spacing, scatterer_grid_size, grid_height);

% Generate A with all the TX/RX antennas
% Alternatively, you can already trim by trimming the vectors
% txPoints/rxPoints 
A = gen_A(lambda, txPoints, rxPoints, scattererPoints);
numRX = size(rxPoints, 1);

% Trim A for selected antennas
txPerm = randperm(rStr, 2 * num_antennas, num_eff_antennas);
rxPerm = randperm(rStr, 2 * num_antennas, num_eff_antennas);
A = trim_A(A, txPerm, rxPerm, numRX);

% Generate inverse 2D Fourier transform matrix 
iF = gen_iF(scatterer_grid_size);

% Get transformed A
AiF = A * iF;
