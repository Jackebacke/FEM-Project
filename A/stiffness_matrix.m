function A=stiffness_matrix(x)
%
% Returns the assembled stiffness matrix A.
% Input is a vector x of node coords.
%

N = length(x) - 1; % number of elements
A = zeros(N+1, N+1); % initialize stiffnes matrix to zero
for i = 1:N % loop over elements
    h = x(i+1) - x(i); % element length
    n = [i i+1]; % nodes
    A(n,n) = A(n,n) + [1 -1; -1 1]/h; % assemble element stiffness
end
% Set dirichlet boundary conditions, by setting first and last row to identity
A(1, 1) = 1;
A(1, 2:end) = 0;
A(N+1, N+1) = 1;
A(N+1, 1:end-1) = 0;
end

