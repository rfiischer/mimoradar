clear;

% Default parameters 
lambda = 0.5;                   % wavelength (cm)
scatterer_grid_size = 32;       % number of scatterers in each axis 
scatterer_spacing = 0.5;        % scatterer spacing (cm)
grid_height = 30;               % grid height (cm)
num_antennas = 46;              % number of antennas in array (on each side)
antenna_spacing = lambda / 2;   % antenna spacing (cm)
sparsity = 0.1;                 % sparsity of scatterer grid (%)
num_eff_antennas = 20;          % effective number of antennas (both sides) chosen randomly 
random_seed = 0;                % 'shuffle' for random; any int for reproducibility

% Random generator
rStr = RandStream('mcg16807', 'Seed', random_seed);

% Get scatterer grid positions 
scattererPoints = xy_grid(scatterer_spacing, ...
    scatterer_spacing, ...
    scatterer_grid_size, ...
    scatterer_grid_size, ...
    grid_height);

% Get TX array positions 
tx1 = xy_grid(0, ...
    antenna_spacing, ...
    1, ...
    num_antennas, ...
    0, ...
    'xCoord', -(num_antennas + 1) / 2 * antenna_spacing);

tx2 = xy_grid(0, ...
    antenna_spacing, ...
    1, ...
    num_antennas, ...
    0, ...
    'xCoord', (num_antennas + 1) / 2 * antenna_spacing);

txPoints = [tx1; tx2];

% Get RX array positions 
rx1 = xy_grid(antenna_spacing, ...
    0, ...
    num_antennas, ...
    1, ...
    0, ...
    'yCoord', -(num_antennas + 1) / 2 * antenna_spacing);

rx2 = xy_grid(antenna_spacing, ...
    0, ...
    num_antennas, ...
    1, ...
    0, ...
    'yCoord', (num_antennas + 1) / 2 * antenna_spacing);

rxPoints = [rx1; rx2];

% Pick random antennas
txPerm = randperm(rStr, 2 * num_antennas, num_eff_antennas);
txPoints = txPoints(txPerm, :);

rxPerm = randperm(rStr, 2 * num_antennas, num_eff_antennas);
rxPoints = rxPoints(rxPerm, :);

% Plot setup
plot_setup(scattererPoints, txPoints, rxPoints);

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

% Compare with function 
[txPointsF, rxPointsF, scattererPointsF] = generate_setup(antenna_spacing, num_antennas, ...
    scatterer_spacing, scatterer_grid_size, grid_height);
Af = gen_A(lambda, txPointsF, rxPointsF, scattererPointsF);
numRX = size(rxPointsF, 1);
Af = trim_A(Af, txPerm, rxPerm, numRX);
fprintf('max|A - Af|: %e\n', max(max(abs(Af - A))));

% Plot correlation of columns
% Here we do PCA in a 3D space (interpreted as color)
B = real_pca(A, 3);
B = color_space(B);
figure;
scatter3(B(:, 1), B(:, 2), B(:, 3));
title('3D Projection of Columns of A');
daspect([1, 1, 1]);
figure;
image(reshape(B, scatterer_grid_size, scatterer_grid_size, []));
title('Color Plot');
daspect([1, 1, 1]);

% Plot correlation in space domain
figure;
ANorm = A ./ vecnorm(A);
pivotI = floor(scatterer_grid_size / 2 + 1/2);
pivotJ = floor(scatterer_grid_size / 2 + 1/2);
pivot = (pivotI - 1) * scatterer_grid_size + pivotJ;
correlations = abs(ANorm(:, pivot)' * ANorm);
imagesc(reshape(correlations, scatterer_grid_size, scatterer_grid_size));
title('Spacial Correlation');
daspect([1, 1, 1]);

% Coherence of A before FFT
fprintf('Coherence of A before ifft: %f\n', coherence(A));

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

% Compare with function
iFf = gen_iF(scatterer_grid_size);
fprintf('max|iF - iFf|: %e\n', max(max(abs(iFf - iF))));

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
AiF = A * iF;

% Plot correlation of columns
% Here we do PCA in a 3D space (interpreted as color)
B = real_pca(AiF, 3);
B = color_space(B);
figure;
scatter3(B(:, 1), B(:, 2), B(:, 3));
title('3D Projection of Columns of A - After IFFT');
daspect([1, 1, 1]);
figure;
image(reshape(B, scatterer_grid_size, scatterer_grid_size, []));
title('Color Plot - After IFFT');
daspect([1, 1, 1]);

% Plot correlation in frequency domain
figure;
ANorm = AiF ./ vecnorm(AiF);
pivotI = floor(scatterer_grid_size / 2 + 1 / 2);
pivotJ = floor(scatterer_grid_size / 2 + 1 / 2);
pivot = (pivotI - 1) * scatterer_grid_size + pivotJ;
correlations = abs(ANorm(:, pivot)' * ANorm);
imagesc(reshape(correlations, scatterer_grid_size, scatterer_grid_size));
title('Frequency Correlation');
daspect([1, 1, 1]);

% Coherence of A after FFT
fprintf('Coherence of A after ifft: %f\n', coherence(AiF));

% Generate signal vector x 
N = size(AiF, 2);
s = floor(sparsity * N);
[x, support] = sparse_x(rStr, N, s, true);

% Sample
y = AiF * x;
[xHat, S] = omp(y, AiF);

% Compare
fprintf('Elements in S^hat but not in S\n')
d1 = setdiff(S, support);
disp(d1)
fprintf('xHat in these entries\n')
disp(xHat(d1))

fprintf('Elements in S but not in S^hat\n')
d2 = setdiff(support, S);
disp(d2)
fprintf('x in these entries\n')
disp(xHat(d2))

fprintf('Error: %e\n', sum(abs(xHat - x) .^ 2))

xHatSpace = ifft2(reshape(xHat, scatterer_grid_size, scatterer_grid_size));
figure
imagesc(abs(xHatSpace))
title('|x| in Space Domain')
daspect([1, 1, 1])
figure
imagesc(angle(xHatSpace))
title('angle(x) in Space Domain')
daspect([1, 1, 1])
figure
imagesc(abs(reshape(xHat, scatterer_grid_size, scatterer_grid_size)))
title('|x| in Frequency Domain')
daspect([1, 1, 1])
