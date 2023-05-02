function [xHat, S] = omp(y, A, tol, rtol)
%OMP Performs OMP given measurements vector y and sensing matrix A
arguments
    y; A;
    tol = 1e-10;
    rtol = 1e-3;
end
% Initialize algorithm
xHat = zeros(size(A, 2), 1);
S = [];

% Begin
i = 0;
while any(abs(y - A * xHat) ./ abs(y) > tol)
    i = i + 1;

    % Get index of smallest error 
    [~, idx] = sort(-abs(A' * (y - A * xHat)));

    % Expand the support
    S = [S; idx(1)];

    % Perform LS approximation
    As = A(:, S(1:i));
    Ast = As' * As;
    if rcond(Ast) < rtol
        break
    end
    At = Ast \ As';
    xHat = zeros(size(xHat));
    xHat(S(1:i)) = At * y;
end
end

