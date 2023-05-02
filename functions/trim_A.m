function A = trim_A(A, txInd, rxInd, numRX)
%TRIM_A Given txInd and rxInd, only select the rows containing these
%antennas

A = A((repelem(txInd, length(rxInd), 1) - 1) * numRX + repmat(rxInd, length(txInd), 1), :);

end

