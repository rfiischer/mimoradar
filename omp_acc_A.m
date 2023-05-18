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
sim_type = 'best';              % choose between 'best' to compute the best
                                % and average performance, or 'random' to
                                % make A random
nSearchIter = 1000;             % how many iterations of random search for good A

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
s = 102;
step = 2;
start = 8;
stop = 32;

if strcmp(sim_type, 'best')
    accuracy = zeros(N, 2);
else
    accuracy = zeros(N, 1);
end

% 'A' database directory
if ~exist('A_database/', 'dir')
    mkdir('A_database/')
end

for Neff = start:step:stop

    right1 = 0;
    right2 = 0;
    M = Neff ^ 2;

    if strcmp(sim_type, 'best')
        fileNameA = fullfile('A_database', sprintf('A_%d.mat', Neff));
        if ~exist(fileNameA, 'file')
            [ACand, AMean] = random_search_A(rStr, AiF, numRX, num_antennas, Neff, nSearchIter, ...
                fileNameA);

        else
            loadA = load(fileNameA);
            ACand = loadA.ACand;
            AMean = loadA.AMean;

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
            [xHat1, ~] = omp(y1, ACand);
            [xHat2, ~] = omp(y2, AMean);
            error1 = sum(abs(xHat1 - x) .^ 2);
            error2 = sum(abs(xHat2 - x) .^ 2);
        
            if error1 < 1e-10
                right1 = right1 + 1;
            end

            if error2 < 1e-10
                right2 = right2 + 1;
            end
        end
    
        accuracy(M, 1) = right1 / nIter;
        accuracy(M, 2) = right2 / nIter;
        save(fileName, 'nIter', 'N', 's', 'accuracy', 'step', 'start', 'stop', 'Neff');

    else

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
                right1 = right1 + 1;
            end
        end
    
        accuracy(M) = right1 / nIter;
        save(fileName, 'nIter', 'N', 's', 'accuracy', 'step', 'start', 'stop');

    end

end
