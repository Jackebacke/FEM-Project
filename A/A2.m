clear; clc; close all;
%% Equation parameters
a = -1; b = 1; % end points of interval
u_a = 0; u_b = 0; % boundary value at x=a and x=b
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


%% Main script

A = stiffness_matrix(x)
B = load_vector(x, @f, u_a, u_b)
M = mass_matrix(x)

hello = 1;