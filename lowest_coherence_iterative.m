clear;

% File parameters
mName = mfilename('fullpath');
fileName = get_name(mName);

% Parameters
M = 3;
N = 4;
nIter = 5000;
alpha = 9e-3;           % gradient scaling factor
tol = 1e-3;             % tolerance on plateau (reduce alpha)
pauseTime = 0.0;        % pause time for animation plot
rStr = RandStream('mcg16807', 'Seed', 0);

% Create A
A = gaussian_A(rStr, M, N, false);
An = A ./ vecnorm(A);

% Plot if 2D or 3D
h = 0;
if M == 2
    h = quiver(zeros(1, N), zeros(1, N), An(1, :), An(2, :));
    xlim([-1, 1]);
    ylim([-1, 1]);
    daspect([1, 1, 1]);
elseif M == 3
    h = quiver3(zeros(1, N), zeros(1, N), zeros(1, N), An(1, :), An(2, :), An(3, :));
    xlim([-1, 1]);
    ylim([-1, 1]);
    zlim([-1, 1]);
    daspect([1, 1, 1]);
end

loss2 = zeros(nIter, 1);
lossMax = zeros(nIter, 1);
for i = 1:nIter
    % Compute correlation matrix
    fprintf('i: %d\n', i);
    C = An' * An;
    C = C .* (ones(size(C)) - eye(size(C)));
    CAbs = abs(C);

    % Compute loss metrics
    loss2(i) = sum(CAbs(:) .^ 2);
    [lossMax(i), idx] = max(CAbs(:));

    % Find columns with largest coherence
    [col1, col2] = ind2sub(size(C), idx);
    signFactor = sign(C(col1, col2));

    % Perform gradient descent
    An(:, col1) = An(:, col1) - alpha * signFactor * An(:, col2);
    An(:, col2) = An(:, col2) - alpha * signFactor * An(:, col1);

    % Normalize A
    An = An ./ vecnorm(An);

    % Reduce alpha if plateau
    if (i > 1000) && (~mod(i, 1000))
        if abs(mean(lossMax(i - 1000:i)) - mean(lossMax(i - 2000 + 1:i - 1000))) < tol
            alpha = alpha / 2;
        end
    end

    % Plot animation if 2D or 3D
    if (M == 2 || M == 3) && ishandle(h) && (pauseTime ~= 0)
        pause(pauseTime);
        h.UData = An(1, :);
        h.VData = An(2, :);
        if N == 3
            h.WData = An(3, :);
        end
    elseif ~ishandle(h)
       break
    end

end

figure;
plot(loss2);

figure;
plot(lossMax);

A = An .* vecnorm(A);
fprintf('Coherence: %.6f\n', coherence(A));

if M == 2
    figure;
    h = quiver(zeros(1, N), zeros(1, N), A(1, :), A(2, :));
    xlim([-1, 1]);
    ylim([-1, 1]);
    daspect([1, 1, 1]);
elseif M == 3
    figure;
    h = quiver3(zeros(1, N), zeros(1, N), zeros(1, N), A(1, :), A(2, :), A(3, :));
    xlim([-1, 1]);
    ylim([-1, 1]);
    zlim([-1, 1]);
    daspect([1, 1, 1]);
end

save(fileName, 'A', 'loss2', 'lossMax');
