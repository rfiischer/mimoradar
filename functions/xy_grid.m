function points = xy_grid(xSpace, ySpace, xNum, yNum, zCoord, xyCoords)
arguments
    xSpace;ySpace;xNum;yNum;zCoord;
    xyCoords.xCoord = 0;
    xyCoords.yCoord = 0;
end
%SCATERER_GRID Creates grid of scaterers
%   The points are then ordered in a column array changing y axis first
xSize = xSpace * (xNum - 1);
ySize = ySpace * (yNum - 1);
xPositions = (0:(xNum - 1))' * xSpace - xSize / 2 + xyCoords.xCoord;
yPositions = (0:(yNum - 1))' * ySpace - ySize / 2 + xyCoords.yCoord;

% For R2023a and above, use 'combinations' function
xyPoints = [repelem(xPositions, yNum, 1), repmat(yPositions, xNum, 1)];
points = [xyPoints, zCoord * ones(xNum * yNum, 1)];

end

