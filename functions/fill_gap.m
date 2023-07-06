function S = fill_gap(S, grid_size)
%FILL_GAP Finds gaps inside sequences in S and fills them

templates = {[2, 2, 1, 2, 2], [2, 2, 1, 2], [2, 1, 2, 2], ...
    [2, 2, 1, 1, 2, 2]};

sPattern = ones(grid_size ^ 2, 1);
sPattern(S) = 2;
sPattern = reshape(sPattern, grid_size, grid_size);

for i = 0:1
    for j = 1:size(templates, 2)
        if i == 0
            template = templates{j}';
        else
            template = templates{j};
        end
        sumTemplate = ones(size(template));
        overhead = find(template == 1, 1) - 1;
        
        normGrid = conv2(sPattern .^ 2, sumTemplate, 'valid');
        patternGrid = conv2(sPattern, flip(template), 'valid') .^ 2 ./ normGrid;
        
        maxVal = sum(template .^ 2);
        newIdx = find(patternGrid(:) == maxVal);
        [row, col] = ind2sub(size(patternGrid), newIdx);
        newIdx = (col + i * overhead - 1) * grid_size + row + (1 - i) * overhead;
        
        S = union(S, newIdx);

    end
end

end

