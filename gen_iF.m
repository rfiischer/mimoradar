function iF = gen_iF(grid_size)
%GEN_IF Generate flattened form of 2D IFFT matrix

% We have that Y = ifft2(X), Y and X are 2D matrices 
% If we flatten X, we can write flatten(Y) = iF * flatten(X)
% Then, each row of iF corresponds to an element Y(p, q), where p changes
% faster
% Analogously, each column of iF corresponds to an element X(j, k), where j
% changes faster

% J is the input X row index, which changes first (y axis)
% K is the input X column index, which changes later (x axis)
% J, K don't change in iF along the rows, just columns
J = repmat((1:grid_size) - 1, grid_size ^ 2, grid_size);
K = repelem((1:grid_size) - 1, grid_size ^ 2, grid_size);

% P, Q don't change in iF along the columns, just rows 
P = J';
Q = K';

exponent = 1i * 2 * pi * (J .* P / grid_size + K .* Q / grid_size);
iF = 1 / grid_size ^ 2 * exp(exponent);

% Alternative way (used to check the above)
% F = zeros(scatterer_grid_size ^ 2, scatterer_grid_size ^ 2);
% for i = 1:scatterer_grid_size
%     for j = 1:scatterer_grid_size
%         x = zeros(scatterer_grid_size, scatterer_grid_size);
%         x(i, j) = 1;
%         F(:, (j-1) * scatterer_grid_size + i) = reshape(fft2(x), [], 1);
%     end
% end
% iF = inv(F);

end

