clear; clc; close all;
%% Equation parameters
a = -1; % left end point of interval
b = 1; % right end point of interval

delta = 0.1;

%% FEM parameters
N = 12;
TOL = 1e-6;
max_nodes = 1e4;
node_start = 12;
lambda = 0.9;
h = (b-a)/N; % initial mesh size
x = a:h:b; % initial node coords

%% Equation functions

function g = g(x)
g = 10 .* x .* sin(7 * pi .* x);
end

function f = f(x)
if g(x) > abs(x)
    f = abs(x);
elseif g(x) < -abs(x)
    f = -abs(x);
else
    f = 10;
end
end

%% FEM functions

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

function B=load_vector(x, u_a, u_b)
%
% Returns the assembled load vector b.
% Input is a vector x of node coords, and boundary values u_a, u_b.
%
N = length(x) - 1;
B = zeros(N+1, 1);
for i = 1:N
    h = x(i+1) - x(i);
    n = [i i+1];
    B(n) = B(n) + [f(x(i)); f(x(i+1))]*h/2; % Trapezoidal rule
end
B(1) = u_a;
B(N+1) = u_b;
end

function M = mass_matrix(x)
%
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
% Set first and last row?

end

%% Main script

A = stiffness_matrix(x)
B = load_vector(x, 0, 0)
M = mass_matrix(x)