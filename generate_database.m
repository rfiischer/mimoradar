clear;

% Default parameters 
lambda = 0.4;                   % wavelength (cm)
scatterer_grid_size = 32;       % number of scatterers in each axis 
scatterer_spacing = 16/31;      % scatterer spacing (cm)
grid_height = 30;               % grid height (cm)
num_antennas = 46;              % number of antennas in array (on each side)
antenna_spacing = 14/45;        % antenna spacing (cm)
random_seed = 'shuffle';        % 'shuffle' for random; any int for reproducibility
nSearchIter = 1000;             % how many iterations of random search for good A
experimentIdx = 4;              % experiment identifier 

start = 8;
step = 2;
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

% 'A' database directory
baseName = sprintf('A_database_%d//', experimentIdx);
if ~exist(baseName, 'dir')
    mkdir(baseName)
end
save(fullfile(baseName, 'parameters.mat'), ...
    'lambda', 'scatterer_grid_size', 'scatterer_spacing', 'grid_height', ...
    'num_antennas', 'antenna_spacing', 'random_seed', 'nSearchIter')

for Neff = start:step:stop
    fprintf('M = %d', Neff ^ 2)

    M = Neff ^ 2;
    fileNameA = fullfile(baseName, sprintf('A_%d.mat', Neff));
    [ACand, AMean] = random_search_A(rStr, AiF, numRX, Neff, nSearchIter, ...
        fileNameA, 0.1, 0.4, 0.4, true, @how_close);

end
