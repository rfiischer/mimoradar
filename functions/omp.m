function [xHat, S] = omp(y, A, tol)
%OMP Performs OMP given measurements vector y and sensing matrix A
arguments
    y; A;
    tol = 1e-10;
end
% Initialize algorithm
xHat = zeros(size(A, 2), 1);
S = [];

% Begin
i = 0;
while any(abs(y - A * xHat) > tol) && (i < size(A, 1))
    i = i + 1;

    % Get index of smallest error 
    [~, idx] = max(abs(A' * (y - A * xHat)));

    % Expand the support
    if ~ismember(idx, S)
        S = [S; idx];
    else
        break
    end

    % Perform LS approximation
    As = A(:, S(1:i));
    
    xHat = zeros(size(xHat));
    xHat(S(1:i)) = pinv(As) * y;
end
end

