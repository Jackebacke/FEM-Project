function M = mass_matrix(x)
% Returns the assembled mass matrix M.
% Input is a vector x of node coords.
%
N = length(x) - 1; % number of elements
M = zeros(N+1, N+1); % initialize stiffnes matrix to zero
for i = 1:N % loop over elements
    h = x(i+1) - x(i); % element length
    n = [i i+1]; % nodes
    M(n,n) = M(n,n) + [1/3 1/6; 1/6 1/3]*h; % assemble element stiffness
end
end
