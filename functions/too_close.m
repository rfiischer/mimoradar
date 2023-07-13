function result = too_close(txPerm, rxPerm, Neff, tol)
%TOO_CLOSE Filter antenna choices that have too many close elements 
arguments
    txPerm; rxPerm; Neff
    tol = 2/5;
end

% Sort permutations 
txPerm = sort(txPerm);
rxPerm = sort(rxPerm);

% Check for consecutive (assuming numerical closeness = spacial closeness)
metric = (sum((txPerm(2:end) - txPerm(1:end-1)) == 1) + ...
    sum((rxPerm(2:end) - rxPerm(1:end-1)) == 1)) / Neff;

% Compare to threshold and decide 
if metric > tol
    result = false;
else
    result = true;
end

end

