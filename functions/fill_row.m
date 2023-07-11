function S = fill_row(S, grid_size, threshold)
%FILL_ROW Finds overfilled rows/columns in OMP and adds the remaining
%entries in the support in hopes LS finds the true entries

sPattern = zeros(grid_size ^ 2, 1);
sPattern(S) = 1;
sPattern = reshape(sPattern, grid_size, grid_size);

rowSum = sum(sPattern, 2);
colSum = sum(sPattern, 1);

rowPer = rowSum(:) / grid_size;
colPer = colSum(:) / grid_size;

populatedRow = find(rowPer > threshold);
populatedCol = find(colPer > threshold);

newRowIdx = repelem(populatedRow, grid_size, 1) + ...
    repmat((0:(grid_size-1))' * grid_size, length(populatedRow), 1);

newColIdx = (repelem(populatedCol, grid_size, 1) - 1) * grid_size + ...
    repmat((1:grid_size)', length(populatedCol), 1);

S = union(S, [newRowIdx; newColIdx]);

end

