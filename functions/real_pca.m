function B = real_pca(A, numComp)
%REAL_PCA Perform PCA of complex vectors to a lower dimensional real vector
%base

% % Both methods are equivalent, the vectors in B may just be -1 multiples of
% % each other
% Ar = real(A);
% Ai = imag(A);
% Ac = [Ar * Ar', Ar * Ai'; Ai * Ar', Ai * Ai'];
% [V, D] = eig(Ac);
% [~, I] = sort(diag(D), 'descend');
% V = V(:, I);
% % B is the PCA projection of the vectors in the columns of A
% B = real(A' * (V(1:size(A, 1), 1:numComp) + 1j * V(size(A, 1) + 1:end, 1:numComp)));

Ar = real(A);
Ai = imag(A);
Ac = Ar' * Ar + Ai' * Ai;
[V, D] = eig(Ac);
[D, I] = sort(diag(D), 'descend');
V = V(:, I);
% B is the PCA projection of the vectors in the columns of A
B = V(:, 1:numComp) .* sqrt(D(1:numComp))';

end

