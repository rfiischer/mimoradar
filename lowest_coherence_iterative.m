% Algorithm to iteratively find a M, N matrix with the smallest coherence
% possible
clear;

% File parameters
mName = mfilename('fullpath');
fileName = get_name(mName);

% Parameters
M = 40;
N = 100;
complex = true;
nIter = 10000;
alpha = 9e-3;           % gradient scaling factor
tol = 1e-3;             % tolerance on plateau (reduce alpha)
pauseTime = 0.0;        % pause time for animation plot
rStr = RandStream('mcg16807', 'Seed', 'shuffle');

% Create A
A = gaussian_A(rStr, M, N, complex);
An = A ./ vecnorm(A);

% Plot if 2D or 3D
h = 0;
if M == 2 && ~complex
    h = quiver(zeros(1, N), zeros(1, N), An(1, :), An(2, :));
    xlim([-1, 1]);
    ylim([-1, 1]);
    daspect([1, 1, 1]);
elseif M == 3 && ~complex
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
    signFactor = angle(C(col1, col2));

    % Perform gradient descent
    An(:, col1) = An(:, col1) - alpha * exp(- 1j * signFactor) * An(:, col2);
    An(:, col2) = An(:, col2) - alpha * exp(+ 1j * signFactor) * An(:, col1);

    % Normalize A
    An = An ./ vecnorm(An);

    % Reduce alpha if plateau
    if (i > 1000) && (~mod(i, 1000))
        if abs(mean(lossMax(i - 1000:i)) - mean(lossMax(i - 2000 + 1:i - 1000))) < tol
            alpha = alpha / 2;
        end
    end

    % Plot animation if 2D or 3D
    if (M == 2 || M == 3) && ~complex && (pauseTime ~= 0)
        pause(pauseTime);
        if ishandle(h)
            h.UData = real(An(1, :));
            h.VData = real(An(2, :));
            if M == 3
                h.WData = real(An(3, :));
            end
        else
            break
        end
    end

end

figure;
plot(loss2);

figure;
plot(lossMax);

A = An .* vecnorm(A);
fprintf('Coherence: %.6f\n', coherence(A));

if M == 2 && ~complex
    figure;
    h = quiver(zeros(1, N), zeros(1, N), real(A(1, :)), real(A(2, :)));
    xlim([-1, 1]);
    ylim([-1, 1]);
    daspect([1, 1, 1]);
elseif M == 3 && ~complex
    figure;
    h = quiver3(zeros(1, N), zeros(1, N), zeros(1, N), real(A(1, :)), real(A(2, :)), real(A(3, :)));
    xlim([-1, 1]);
    ylim([-1, 1]);
    zlim([-1, 1]);
    daspect([1, 1, 1]);
end

save(fileName, 'A', 'loss2', 'lossMax');
