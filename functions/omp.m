function [xHat, S] = omp(y, A, post_process, tol, rtol)
%OMP Performs OMP given measurements vector y and sensing matrix A
arguments
    y; A;
    post_process = @(x, y)(x);
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
    [zHat, stop] = ls_approximation(A, S, y, rtol);
    if stop
        S = S(1:end-1);
        break;
    else
        xHat = zeros(size(A, 2), 1);
        xHat(S) = zHat;
    end

end

% Apply post processing into the estimated support 
S = post_process(S, sqrt(size(A, 2)));
[zHat, ~] = ls_approximation(A, S, y, 0);
xHat = zeros(size(A, 2), 1);
xHat(S) = zHat;

end

