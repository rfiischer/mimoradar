function S = fill_gap(S, grid_size)
%FILL_GAP Finds gaps inside sequences in S and fills them

templates = {[2, 2, 2, 1], [1, 2, 2, 2], [2, 2, 1, 2, 2], ...
    [2, 2, 1, 2], [2, 1, 2, 2], [2, 2, 1, 1, 2, 2]};

% Lay out the support in a square grid shape
sPattern = ones(grid_size ^ 2, 1);
sPattern(S) = 2;
sPattern = reshape(sPattern, grid_size, grid_size);

% i == 0 is vertical
% i == 1 is horiziontal 
for i = 0:1
    for j = 1:size(templates, 2)
        if i == 0
            template = templates{j}';
        else
            template = templates{j};
        end
        sumTemplate = ones(size(template));
        gaps = find(template == 1);             % index of gap 
        
        % This sliding window gives the norm squared of each sequence to be
        % compared 
        normGrid = conv2(sPattern .^ 2, sumTemplate, 'valid');
        patternGrid = conv2(sPattern, flip(template), 'valid') .^ 2 ./ normGrid;
        
        maxVal = sum(template .^ 2);
        newIdx = find(patternGrid(:) == maxVal);
        [row, col] = ind2sub(size(patternGrid), newIdx);
        numToAdd = size(newIdx, 1);

        % Fill in the gap
        gapsize = length(gaps);
        newIdxLinear = zeros(1, numToAdd * gapsize);
        for k = 1:gapsize
            % Convert row, col to a linear index 
            overhead = gaps(k) - 1;

            % If i == 0, then just the row index needs correction due to
            % the 'valid' convolution
            % If i == 1, then just the column index needs correction
            newIdxLinear((k - 1) * numToAdd + 1:k * numToAdd) = (col - 1 + i * overhead) * grid_size + ...
                row + (1 - i) * overhead;
        end
        
        % Update S and regenerate pattern 
        S = union(S, newIdxLinear);
        sPattern = ones(grid_size ^ 2, 1);
        sPattern(S) = 2;
        sPattern = reshape(sPattern, grid_size, grid_size);

    end
end

end

