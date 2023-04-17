clear;

% Random generator
rStr = RandStream('mcg16807', 'Seed', 0);

nIter = 1000;
N = 1000;
M = 300;
s = 10;

right = 0;

for i = 1:nIter
    % Generate sparse x
    x = zeros(N, 1);
    support = randperm(rStr, N, s);
    x(support) = randn(rStr, s, 1);

    % Generate A
    A = 1 / sqrt(M) * randn(rStr, M, N);
    y = A * x;

    [~, S] = omp(y, A);

    if length(S) == length(support)
        if all(sort(find(x)) == sort(S))
            right = right + 1;
        end
    end
end

% Coherence of A
fprintf('Coherence of A: %f\n', coherence(A));
