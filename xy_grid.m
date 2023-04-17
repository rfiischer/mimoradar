function points = xy_grid(xSize, ySize, xNum, yNum, zCoord, xyCoords)
arguments
    xSize;ySize;xNum;yNum;zCoord;
    xyCoords.xCoord = 0;
    xyCoords.yCoord = 0;
end
%SCATERER_GRID Creates grid of scaterers
%   The points are then ordered in a column array changing y axis first
xPositions = linspace(0, xSize, xNum)' - xSize / 2 + xyCoords.xCoord;
yPositions = linspace(0, ySize, yNum)' - ySize / 2 + xyCoords.yCoord;

% For R2023a and above, use 'combinations' function
xyPoints = [repelem(xPositions, yNum, 1), repmat(yPositions, xNum, 1)];
points = [xyPoints, zCoord * ones(xNum * yNum, 1)];

end

