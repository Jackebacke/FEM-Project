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
xnodes= a:h:b; % initial node coords

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

function val = phi_i(x, i, xnodes)
% phi_i(x, i, xnodes)
% Returns value of the i-th hat function at point(s) x
%
% xnodes : vector of mesh nodes [x0, x1, ..., xN]
h = xnodes(2) - xnodes(1);
    if i == 1
        val = x/h;
    elseif i == numel(xnodes)
        val = (xnodes(i-1)-x)/h;
    else
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

end

%% Skapa vektor av basfunktionerna
function Phi = basvec(x, xnodes)
    % Define nodes
% Initialize matrix to store values of each phi
N = numel(xnodes);
Phi = zeros(numel(xnodes),1);
% Compute phi_i at x
for i = 1:N               % exclude first and last if you use boundary conditions
    Phi(i) = phi_i(x, i, xnodes);
end
end

function Laplu_h = Laplu_h(x, basvec,Z)
    disp(size(Z))
    disp(size(basvec))
    Laplu_h = Z'* basvec;
end
%% Main script

% while sum(eta2) > TOL && length(x) < max_nodes

% end

A = stiffness_matrix(xnodes);
B = load_vector(xnodes, @f, u_a, u_b);
M = mass_matrix(xnodes);

u_h = fem_solver(xnodes, u_a, u_b);
Z = -M \ A*u_h;




% function Phi = basvec(x, xnodes)
%     % Define nodes
% N = 10;                     % number of nodes
% xnodes = linspace(0, 1, N); % nodes from 0 to 1

% % Define evaluation points (can be finer than nodes)
% x = linspace(0, 1, 100);    % 100 points for plotting

% % Initialize matrix to store values of each phi_i
% Phi = zeros(numel(x), numel(xnodes));

% % Compute phi_i at all x
% for i = 2:N-1               % exclude first and last if you use boundary conditions
%     Phi(:, i) = phi_i(x, i, xnodes);
% end
% end
% basvec =
% x = linspace(0, 0.9, 10);

phi = basvec(0.1, xnodes);

disp(Laplu_h(0.1,phi,Z))


%%% INTEGRTION FÖR FEL APPROXIMERINGEN

function trapezoidal = num_integration(xnodes)
n = numel(xnodes);
h = xnodes(2) -xnodes(1);
y_values = zeros(numel(xnodes), 1);
for i = 1:numel(xnodes)             
    y_values(i) = phi_i(xnodes(i), i, xnodes) + f(xnodes(i));
end
trapezoidal = h*(1/2*(y_values(1)+ y_values(n)) + sum(y_values(2:n-1)));
end
disp(num_integration(xnodes))


eta_i = zeros(N,1);
for i = 1:N
    eta_i(i) = num_integration([xnodes(i) xnodes(i+1)]);
end
disp(eta_i)
disp(size(eta_i))

%%FEL APPROXIMERING
fel = h*sum(eta_i);
disp(fel)