function [xHat, stop] = ls_approximation(A, S, y, rtol)
%LS_APPROXIMATION Gets the least squares solution of an overdetermined
%system of equations

As = A(:, S);

Ast = As' * As;
if rcond(Ast) < rtol
    stop = true;
    xHat = 0;
else
    stop = false;
    At = Ast \ As';
    xHat = At * y;
end

end

