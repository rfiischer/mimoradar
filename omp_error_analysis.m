clear;

% Default parameters 
random_seed = 'shuffle';        % 'shuffle' for random; any int for reproducibility
nIter = 1000;                   % number of attempts 
directory = 'A_database//';     % experiment identifier 
Neff = 20;                      % effective number of antennas
s = 0.1;                        % sparsity
grid_size = 32;                 % scatterer grid size
animation_start_idx = 0;
post_process = @fill_gap;       % post processing function

% Random generator
rStr = RandStream('mcg16807', 'Seed', random_seed);

% Load A
M = Neff ^ 2;
fileNameA = fullfile(directory, sprintf('A_%d.mat', Neff));
loadA = load(fileNameA);
A = loadA.ACand;

N = size(A, 2);
for i = 1:nIter
    % Generate sparse x
    [x, supp] = sparse_x(rStr, N, s, true);
    
    % Sample 
    y = A * x;
    
    % Perform OMP
    [xHat, suppHat] = omp(y, A, post_process);
    error = sum(abs(xHat - x) .^ 2);
    
    fprintf('Error = %f\n', error);

    suppOriginal = suppHat;
    f1 = figure;
    f2 = figure;

    for j = (length(suppOriginal) + animation_start_idx):length(suppOriginal)
        suppHat = suppOriginal(1:j);

        figure(f1);
        ax = axes;
        [X, Y] = meshgrid(1:grid_size);
        imagesc(ax, reshape(abs(xHat - x), grid_size, grid_size));
        colormap(f1, 'hot');
        hold on;
    
        leftOut = supp(~ismember(supp, suppHat));
        right = supp(ismember(supp, suppHat));
        extra = suppHat(~ismember(suppHat, supp));
    
        sc = scatter(X(leftOut), Y(leftOut), 100, 's', 'DisplayName', '$S \notin \hat{S}$');
        sc.MarkerEdgeColor = [1, 0, 0];
        sc.LineWidth = 2;
    
        sc = scatter(X(right), Y(right), 100, 's', 'DisplayName', '$S \in \hat{S}$');
        sc.MarkerEdgeColor = [0, 1, 0];
        sc.LineWidth = 2;
    
        sc = scatter(X(extra), Y(extra), 100, 's', 'DisplayName', '$\hat{S} \notin S$');
        sc.MarkerEdgeColor = [0, 0, 1];
        sc.LineWidth = 2;
    
        daspect([1, 1, 1]);
        h = legend('show');
        set(h, 'Interpreter','latex');
        title('$|\hat{x} - x|$', 'Interpreter', 'latex');
        set(ax, 'YDir', 'normal');
        colorbar;
    
        if leftOut
            figure(f2);
            ax = axes;
            ANorm = A ./ vecnorm(A);
            C = log(abs(reshape(ANorm(:, leftOut(1))' * ANorm, grid_size, grid_size)));
            imagesc(ax, C);
            title('$\log{|a_e^H A|}$', 'Interpreter', 'latex');
            hold on;
        
            sc = scatter(ax, X(leftOut), Y(leftOut), 100, 's', 'DisplayName', '$S \notin \hat{S}$');
            sc.MarkerEdgeColor = [1, 0, 0];
            sc.LineWidth = 2;
        
            sc = scatter(ax, X(right), Y(right), 100, 's', 'DisplayName', '$S \in \hat{S}$');
            sc.MarkerEdgeColor = [0, 1, 0];
            sc.LineWidth = 2;
        
            sc = scatter(ax, X(extra), Y(extra), 100, 's', 'DisplayName', '$\hat{S} \notin S$');
            sc.MarkerEdgeColor = [0, 0, 1];
            sc.LineWidth = 2;
        
            daspect([1, 1, 1]);
            h = legend('show');
            set(h, 'Interpreter', 'latex');
            set(ax, 'YDir', 'normal');
        end
    
        input('');
    
    end

    waitfor(f1);
    waitfor(f2);

end
