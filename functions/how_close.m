function metric = how_close(txPerm)
%HOW_CLOSE Compute how many adjacent elements there are in array 

% Sort permutations 
txPerm = sort(txPerm);

% Check for consecutive (assuming numerical closeness = spacial closeness)
metric = sum((txPerm(2:end) - txPerm(1:end-1)) == 1);

end

