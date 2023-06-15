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

% We use 'all' instead of 'any' because usually when OMP finds the right
% solution, in the same iteration that all(abs(y - A * xHat) > tol) becomes
% false so does the statement any(abs(y - A * xHat) > tol), that means as
% soon as one entry is minimized so are the others

% Doing this we can stop OMP at fewer iterations, since for the cases a
% mistake happens it is not worth checking until any(...) is false

while all(abs(y - A * xHat) > tol) && (i < size(A, 1))
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
    At = Ast \ As';
    xHat(S(1:i)) = At * y;

end
end

