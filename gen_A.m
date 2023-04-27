function A = gen_A(lambda, size_of_hand, scatterer_grid_size, grid_height, num_antennas)
%GEN_A Generate sampling matrix A

% Get scatterer grid positions 
scattererPoints = xy_grid(size_of_hand, ...
    size_of_hand, ...
    scatterer_grid_size, ...
    scatterer_grid_size, ...
    grid_height);

% Get TX array positions 
tx1 = xy_grid(0, ...
    (num_antennas - 1) * lambda / 2, ...
    1, ...
    num_antennas, ...
    0, ...
    "xCoord", -(num_antennas + 1) / 2 * lambda / 2);

tx2 = xy_grid(0, ...
    (num_antennas - 1) * lambda / 2, ...
    1, ...
    num_antennas, ...
    0, ...
    "xCoord", (num_antennas + 1) / 2 * lambda / 2);

txPoints = [tx1; tx2];

% Get RX array positions 
rx1 = xy_grid((num_antennas - 1) * lambda / 2, ...
    0, ...
    num_antennas, ...
    1, ...
    0, ...
    "yCoord", -(num_antennas + 1) / 2 * lambda / 2);

rx2 = xy_grid((num_antennas - 1) * lambda / 2, ...
    0, ...
    num_antennas, ...
    1, ...
    0, ...
    "yCoord", (num_antennas + 1) / 2 * lambda / 2);

rxPoints = [rx1; rx2];

% Get distances
numScatterers = size(scattererPoints, 1);
numTX = size(txPoints, 1);
numRX = size(rxPoints, 1);

vecTXtoScat = repmat(txPoints, numScatterers, 1) - repelem(scattererPoints, numTX, 1);
dTXtoScat = vecnorm(vecTXtoScat, 2, 2);
dTXtoScat = reshape(dTXtoScat, numTX, []);

vecRXtoScat = repmat(rxPoints, numScatterers, 1) - repelem(scattererPoints, numRX, 1);
dRXtoScat = vecnorm(vecRXtoScat, 2, 2);
dRXtoScat = reshape(dRXtoScat, numRX, []);

% Get vectors a and b
% Each column represents one scatterer 
% Each row represents one TX/RX 
a = exp(1i * 2 * pi * dTXtoScat / lambda) ./ (dTXtoScat .^ 2);
b = exp(1i * 2 * pi * dRXtoScat / lambda) ./ (dRXtoScat .^ 2);

% Assemble measurement matrix
A = zeros(numTX * numRX, numScatterers);
for i = 1:size(scattererPoints, 1)
    A(:, i) = kron(a(:, i), b(:, i));

end

% Each column is one TX
A = reshape(A, numTX, numRX, []);

end
