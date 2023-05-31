clear;

% Random generator
rStr = RandStream('mcg16807', 'Seed', 0);

nIter = 1000;
N = 1024;
M = 400;
s = 0.1;

right = 0;

for i = 1:nIter
    % Generate sparse x
    [x, support] = sparse_x(rStr, N, s, true);

    % Generate A
    A = gaussian_A(rStr, M, N, true);
    y = A * x;

    % Perform OMP
    [xHat, S] = omp(y, A);
    error = sum(abs(xHat - x) .^ 2);

    if all(ismember(support, S))
        if error < 1e-10
            right = right + 1;
        end
    end
end

fprintf('Coherence of A: %f\n', coherence(A));
fprintf('Percentage of right: %2.2f\n', 100 * right / nIter);

