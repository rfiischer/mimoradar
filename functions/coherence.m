function c = coherence(A)
%COHERENCE Returns the coherence of matrix A

% Find norm of columns
normColumns = sqrt(sum(abs(A) .^ 2));

% Normalize columns
A = A ./ normColumns;

% Contains in element (i, j) the product ai' * aj
coherenceMatrix = A' * A;

coherenceCandidates = abs(coherenceMatrix .* (1 - eye(size(coherenceMatrix, 1))));
c = max(coherenceCandidates(:));
end

