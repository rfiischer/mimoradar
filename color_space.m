function B = color_space(B)
%COLOR_SPACE Normalize 3D vectors to fit the color cube
B = B - min(B);
B = B ./ max(B);
end

