clear;

% Default parameters 
lambda = 0.5;                   % wavelength (cm)
size_of_hand = 16;              % hand box size (cm)
scatterer_grid_size = 32;       % number of scatterers in each axis 
grid_height = 100;              % grid height (cm)
num_antennas = 46;              % number of antennas in array (on each side)
sparsity = 0.1;                 % sparsity of scatterer grid (%)
num_eff_antennas = 20;          % effective number of antennas (both sides) chosen randomly 

% Random generator
rStr = RandStream('mcg16807', 'Seed', 0);

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

% Pick random antennas
txPerm = randperm(rStr, 2 * num_antennas, num_eff_antennas);
txPoints = txPoints(txPerm, :);

rxPerm = randperm(rStr, 2 * num_antennas, num_eff_antennas);
rxPoints = rxPoints(rxPerm, :);

% Plot setup
hold on;
scatter3(scattererPoints(:, 1), scattererPoints(:, 2), scattererPoints(:, 3), "blue", "*", "DisplayName", "Scatterers");
scatter3(txPoints(:, 1), txPoints(:, 2), txPoints(:, 3), "green", "s", "DisplayName", "TX");
scatter3(rxPoints(:, 1), rxPoints(:, 2), rxPoints(:, 3), "red", "x", "DisplayName", "RX");
hold off;
view(45, 45);
legend;
xlabel('x (cm)');
ylabel('y (cm)');
zlabel('z (cm)');
daspect([1, 1, 1]);

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

% We have that Y = ifft2(X), Y and X are 2D matrices 
% If we flatten X, we can write flatten(Y) = iF * flatten(X)
% Then, each row of iF corresponds to an element Y(p, q), where p changes
% faster
% Analogously, each column of iF corresponds to an element X(j, k), where j
% changes faster

% J is the input X row index, which changes first (y axis)
% K is the input X column index, which changes later (x axis)
% J, K don't change in iF along the rows, just columns
J = repmat((1:scatterer_grid_size) - 1, scatterer_grid_size ^ 2, scatterer_grid_size);
K = repelem((1:scatterer_grid_size) - 1, scatterer_grid_size ^ 2, scatterer_grid_size);

% P, Q don't change in iF along the columns, just rows 
P = J';
Q = K';

exponent = 1i * 2 * pi * (J .* P / scatterer_grid_size + K .* Q / scatterer_grid_size);
iF = 1 / scatterer_grid_size ^ 2 * exp(exponent);

% Alternative way (used to check the above)
% F = zeros(scatterer_grid_size ^ 2, scatterer_grid_size ^ 2);
% for i = 1:scatterer_grid_size
%     for j = 1:scatterer_grid_size
%         x = zeros(scatterer_grid_size, scatterer_grid_size);
%         x(i, j) = 1;
%         F(:, (j-1) * scatterer_grid_size + i) = reshape(fft2(x), [], 1);
%     end
% end
% iF2 = inv(F);

% Get sampling matrix (after ifft)
A = A * iF;
