function [xHat, S] = omp(y, A, tol, rtol)
%OMP Performs OMP given measurements vector y and sensing matrix A
arguments
    y; A;
    tol = 1e-5;
    rtol = 1e-3;
end
% Initialize algorithm
xHat = zeros(size(A, 2), 1);
S = [];

% Begin
i = 0;
while any(abs(y - A * xHat) > tol)
    i = i + 1;

    % Get index of smallest error 
    [~, idx] = max(abs(A' * (y - A * xHat)));

    % Expand the support
    if ~ismember(idx, S)
        S = [S; idx];
    else
        % This case happens when all the elements of abs(A' * (y - A *
        % xHat)) are almost zero, so it has essentially found the LS
        % solution and 'tol' is too small to be met 
        break
    end

    % Perform LS approximation
    As = A(:, S(1:i));
    xHat = zeros(size(xHat));
    
    Ast = As' * As;
    if rcond(Ast) < rtol
        break
    end

    At = Ast \ As';
    xHat(S(1:i)) = At * y;

end
end

