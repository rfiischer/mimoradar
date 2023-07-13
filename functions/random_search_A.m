function [ACand, AMean, co, minCand, meanCand, txPerm, rxPerm, i] = random_search_A(rStr, AiF, ...
    numRX, num_antennas, num_eff_antennas, numIter, fileName, verbose, processing_function, ...
    tol)
%RANDOM_SEARCH_A Find A with good coherence based on random search
arguments
    rStr; AiF; numRX; num_antennas; num_eff_antennas; numIter;
    fileName = false;
    verbose = false;
    processing_function = @(x, y, z)(true);
    tol = 0.005;
end

co = zeros(numIter, 1);
minCand = 1;
meanCand = 0;
ACand = 0;
AMean = 0;
meanCo = 0;
for i = 1:numIter
    if verbose
        fprintf('i: %d\n', i);
    end

    % Pick random antennas
    txPerm = randperm(rStr, 2 * num_antennas, num_eff_antennas);
    rxPerm = randperm(rStr, 2 * num_antennas, num_eff_antennas);

    % Get subset of A
    ATrim = trim_A(AiF, txPerm, rxPerm, numRX);

    % Compute coherence
    co(i) = coherence(ATrim);

    % Check if it is near mean
    if abs(co(i) - meanCo) < tol
        AMean = ATrim;
        meanCand = co(i);
    end

    % Update mean 
    meanCo = (meanCo * (i - 1) + co(i)) / i;
    
    % Check for the best 
    if (co(i) < minCand) && processing_function(txPerm, rxPerm, num_eff_antennas)
        minCand = co(i);
        ACand = ATrim;
        if fileName
            save(fileName, 'num_antennas', 'num_eff_antennas', 'numIter', 'co', 'minCand', 'ACand', ...
                'txPerm', 'rxPerm', 'i', 'meanCand', 'AMean');
        end
    end
    
end
end

