function A = gaussian_A(rStr, M, N, complex)
arguments
    rStr; M; N;
    complex = false;
end
%GAUSSIAN_A Generate Gaussian sampling matrix 
%   Detailed explanation goes here
if complex
    A = 1 / sqrt(2 * M) * (randn(rStr, M, N) + 1j * randn(rStr, M, N));
else
    A = 1 / sqrt(M) * randn(rStr, M, N);
end
end

