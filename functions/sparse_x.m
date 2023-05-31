function [x, support] = sparse_x(rStr, N, s, complex)
%SPARSE_X Generate sparse vector x using Bernoulli-Gaussian model 
arguments
    rStr; N; s;
    complex = false;
end
x = zeros(N, 1);
support = find(rand(rStr, N, 1) < s);
suppSize = size(support, 1);
if complex
    x(support) = 1 / sqrt(2) * (randn(rStr, suppSize, 1) + 1i * randn(rStr, suppSize, 1));
else
    x(support) = randn(rStr, suppSize, 1);
end
end

