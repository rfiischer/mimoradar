function A = trim_A(A, txInd, rxInd)
%TRIM_A Given txInd and rxInd, only select the rows containing these
%antennas

A = A(rxInd, txInd, :);
A = reshape(A, [], size(A, 3));

end

