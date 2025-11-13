clear; clc; close all;
%% Equation parameters
a = -1; b = 1; % end points of interval
u_a = 0; u_b = 0; % boundary value at x=a and x=b
delta = 0.1;

%% FEM parameters
N = 12;
TOL = 1e-3;
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
    f = g(x);
end
end

%% FEM functions

function u_h = fem_solver(x, u_a, u_b)
A = stiffness_matrix(x);
B = load_vector(x, @f, u_a, u_b);
u_h = A\B; % solve system of equations
end


%% Main script

% while sum(eta2) > TOL && length(x) < max_nodes

% end

A = stiffness_matrix(x)
B = load_vector(x, @f, u_a, u_b)
M = mass_matrix(x)

u_h = fem_solver(x, u_a, u_b)
Z = -M \ A*u_h


%% Basfunktioner  (phi)
function val = phi_i(x, i, xnodes)
% phi_i(x, i, xnodes)
% Returns value of the i-th hat function at point(s) x
%
% xnodes : vector of mesh nodes [x0, x1, ..., xN]

    xi_minus = xnodes(i-1);
    xi       = xnodes(i);
    xi_plus  = xnodes(i+1);
    h = xnodes(2) - xnodes(1);  % assumes uniform spacing

    val = zeros(size(x));  % initialize output

    % Left linear segment
    left_mask = (x >= xi_minus) & (x <= xi);
    val(left_mask) = (x(left_mask) - xi_minus) / h;

    % Right linear segment
    right_mask = (x >= xi) & (x <= xi_plus);
    val(right_mask) = (xi_plus - x(right_mask)) / h;
end
%% Skapa vektor av basfunktionerna
function Phi = basvec(x, xnodes)
    % Define nodes
N = 10;                     % number of nodes
xnodes = linspace(0, 1, N); % nodes from 0 to 1

% Define evaluation points (can be finer than nodes)
x = linspace(0, 1, 100);    % 100 points for plotting

% Initialize matrix to store values of each phi_i
Phi = zeros(length(x), numel(xnodes));

% Compute phi_i at all x
for i = 2:N-1               % exclude first and last if you use boundary conditions
    Phi(:, i) = phi_i(x, i, xnodes);
end
end
x = linspace(0, 0.9, 10);


Laplu_h = 
function I = trapz_integration(f, a, b, n)
%TRAPZ_INTEGRATION  Numerical integration using the trapezoidal rule.
%
%   I = trapz_integration(f, a, b, n)
%
%   Inputs:
%     f - function handle (e.g. @(x) x.^2)
%     a - lower limit of integration
%     b - upper limit of integration
%     n - number of subintervals
%
%   Output:
%     I - approximate value of the integral ∫_a^b f(x) dx
% Define g(x)
g = @(x) 10 .* x .* sin(7 * pi .* x);

% Define f(x) using elementwise and logical operations
f = @(x) (g(x) > abs(x)) .* abs(x) + ...
         (g(x) < -abs(x)) .* (-abs(x)) + ...
         (g(x) <= abs(x) & g(x) >= -abs(x)) .* g(x);
h = (b - a) / n;            % step size
x = a:h:b;                  % grid points
y = f(x) + ;                   % function values
I = (h/2) * (y(1) + 2*sum(y(2:end-1)) + y(end));  % trapezoidal formula
end

eta2 = zeros(N,1);
for i = 1:N
    h = x(i+1) - x(i);
    eta2(i) = 1;
end

