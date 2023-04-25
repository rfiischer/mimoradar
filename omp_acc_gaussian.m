clear;

% Random generator
rStr = RandStream('mcg16807', 'Seed', 0);

nIter = 1000;
N = 1024;
s = 102;

accuracy = zeros(N, 1);
for M = 1:N

    right = 0;

    for i = 1:nIter
        % Progress
        fprintf('Iter: %d, M = %d\n', i, M);

        % Generate sparse x
        [x, support] = sparse_x(rStr, N, s, true);
    
        % Generate A
        A = gaussian_A(rStr, M, N, true);
        y = A * x;
    
        % Perform OMP
        [xHat, S] = omp(y, A);
        error = sum(abs(xHat - x) .^ 2);
    
        if error < 1e-10
            right = right + 1;
        end
    end

    accuracy(M) = right / nIter;

end

if ~exist('results', 'dir')
   mkdir('results');
end

mName = mfilename('fullpath');
fileName = get_name(mName);
save(fileName, 'nIter', 'N', 's', 'accuracy');
