function x = fixedPointIteration(F, x0, tol, maxIter)
% fixedPointIteration - Solves x = F(x) using fixed-point iteration
% Usage: x = fixedPointIteration(F, x0, tol, maxIter)
% Inputs:
%   F       - Function handle for the iteration function F(x)
%   x0      - Initial guess (column vector)
%   tol     - Tolerance for convergence
%   maxIter - Maximum number of iterations
% Output:
%   x       - Approximate solution

x = x0;
for iter = 1:maxIter
    x_new = F(x);
    if norm(x_new - x, 2) < tol
        fprintf('Converged in %d iterations.\n', iter);
        x = x_new;
        return;
    end
    x = x_new;
end

warning('Maximum iterations reached without convergence.');
end


