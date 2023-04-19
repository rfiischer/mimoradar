function [x, support] = sparse_x(rStr, N, s, complex)
%SPARSE_X Generate sparse vector x using Bernoulli-Gaussian model 
arguments
    rStr; N; s;
    complex = false;
end
x = zeros(N, 1);
support = randperm(rStr, N, s);
if complex
    x(support) = 1 / sqrt(2) * (randn(rStr, s, 1) + 1i * randn(rStr, s, 1));
else
    x(support) = randn(rStr, s, 1);
end
end

