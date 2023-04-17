function [xHat, S] = omp(y, A)
%OMP Performs OMP given measurements vector y and sensing matrix A

% Initialize algorithm
xHat = zeros(size(A, 2), 1);
S = [];

% Begin
i = 0;
while y ~= A * xHat
% while sum(abs(y - A * xHat) .^ 2) > 1e-10
    i = i + 1;

    % Get index of smallest error 
    [~, idx] = sort(-abs(A' * (y - A * xHat)));

    % Expand the support
    S = [S; idx(1)];

    % Perform LS approximation
    As = A(:, S(1:i));
    At = (As' * As) \ As';
    xHat = zeros(size(xHat));
    xHat(S(1:i)) = At * y;
end
end

