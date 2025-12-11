clear; clc; close all;
%% Params
T = 2;
dt = 0.01;

d1 = 0.01;
alpha = 4;

function v0 = initial_condition(x,y)
w = rand(1,1);
v0 = 1 + 20 * w;
end

function s = N(u, alpha)
s = -u .* (1 - u) + u ./ (alpha + u);
end

%% Geometry
geometry = @circleg ;
hmax = 1;
[p , e , t] = initmesh(geometry ,'hmax' , hmax);
[A,M,b] = assema(p,t, 1, 1, 0);
u0 = initial_condition(p(1,:), p(2,:))';

%% Time-stepping using fixed-point iteration
psi_current = u0;
tol = 1e-3;
maxIter = 10000;

% Create wrapper function that captures all needed variables
F_wrapped = @(psi_next) M * ((psi_next - psi_current)/dt + N((psi_next + psi_current)/2, alpha)) - d1 * A * ((psi_next + psi_current)/2) + psi_next;

psi_next = fixedPointIteration(F_wrapped, psi_current, tol, maxIter)






