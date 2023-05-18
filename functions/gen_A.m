function A = gen_A(lambda, ...
    txPoints, rxPoints, scattererPoints)
%GEN_A Generate sampling matrix A

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

end

