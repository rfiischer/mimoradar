function [ACand, AMean] = random_search_A(rStr, AiF, ...
    numRX, num_eff_antennas, numIter, fileName, ...
    q1, q2, q3, verbose, metric_function)
%RANDOM_SEARCH_A Find A with good coherence based on random search
arguments
    rStr; AiF; numRX; num_eff_antennas; numIter;
    fileName = false;
    q1 = 0; q2 = 0; q3 = 0;
    verbose = false;
    metric_function = @(x)(0);
end

co = zeros(numIter, 1);
metric = zeros(numIter, 2);
minCoCand = 0;
minM1Cand = 0;
minM2Cand = 0;
meanCoCand = 0;
meanM1Cand = 0;
meanM2Cand = 0;
meanCo = 0;
meanM1 = 0;
meanM2 = 0;
ACand = 0;
AMean = 0;
for i = 1:numIter
    if verbose
        fprintf('i: %d\n', i);
    end

    % Pick random antennas
    txPerm = randperm(rStr, numRX, num_eff_antennas);
    rxPerm = randperm(rStr, numRX, num_eff_antennas);

    % Get subset of A
    ATrim = trim_A(AiF, txPerm, rxPerm, numRX);

    % Compute coherence
    coValue = coherence(ATrim);
    co(i) = coValue;

    % Compute metric
    m1Value = metric_function(txPerm);
    m2Value = metric_function(rxPerm);
    metric(i, 1) = m1Value;
    metric(i, 2) = m2Value;

    % Check if it is near mean
    [~, idx] = min(sqrt((co(1:i) - meanCo) .^ 2 + ...
        (metric(1:i, 1) - meanM1) .^ 2 + ...
        (metric(1:i, 2) - meanM2) .^ 2));

    if idx == i
        AMean = ATrim;
        meanCoCand = coValue;
        meanM1Cand = m1Value;
        meanM2Cand = m2Value;
    end

    % Update mean 
    meanCo = (meanCo * (i - 1) + coValue) / i;
    meanM1 = (meanM1 * (i - 1) + m1Value) / i;
    meanM2 = (meanM2 * (i - 1) + m2Value) / i;
    
    % Check for coherence and metrics quantiles 
    coQuantile = quantile(co(1:i), q1);
    m1Quantile = quantile(metric(1:i, 1), q2);
    m2Quantile = quantile(metric(1:i, 2), q3);

    % Check for the best
    
    if (coValue <= coQuantile) && (m1Value <= m1Quantile) && (m2Value <= m2Quantile)
        minCoCand = coValue;
        minM1Cand = m1Value;
        minM2Cand = m2Value;
        ACand = ATrim;
    end
    
end

if fileName
    save(fileName, 'numRX', 'num_eff_antennas', 'numIter', ...
        'q1', 'q2', 'q3', 'metric_function', ...
        'minCoCand', 'minM1Cand', 'minM2Cand', ...
        'meanCoCand', 'meanM1Cand', 'meanM2Cand', ...
        'ACand', 'AMean', 'i');
end

end

