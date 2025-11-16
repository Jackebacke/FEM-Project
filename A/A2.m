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
function y = f(x_vals)
gx = 10 .* x_vals .* sin(7 * pi .* x_vals);

y = zeros(size(x_vals));
for idx = 1:length(x_vals)
    x = x_vals(idx);
    g = gx(idx);

    if g > abs(x)
        y(idx) = abs(x);
    elseif g < -abs(x)
        y(idx) = -abs(x);
    else
        y(idx) = g;
    end
end
end

%% FEM functions

%% Main script
eta2 = ones(N,1); % Initialize error estimator

% while sum(eta2) > TOL && length(x) < max_nodes
A = stiffness_matrix(x);
B = load_vector(x, @f, u_a, u_b);
u_h = A\B; % solution

M = mass_matrix(x)
Zeta = -M \ A*u_h % Second derivative approximation at nodes

eta2 = zeros(N,1);
for i = 1:N
    h = x(i+1) - x(i);
    eta2(i) = h *1; %TODO: implement error estimator
end
% end


%% Tror inte reesten behövs?

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



Lapl_u_h = 3
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
y = f(x) + 1;                   % function values
I = (h/2) * (y(1) + 2*sum(y(2:end-1)) + y(end));  % trapezoidal formula
end


