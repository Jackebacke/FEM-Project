clear; clc; close all;
%% Params
T = 2;
dt = 0.01;
d1 = 0.1;
alpha = 10;

tol = 1e-3;
maxIter = 1000;

tvec = 0:dt:T

function v0 = initial_condition(x, y)
% Generate random initial condition for each point
w = rand(size(x));
v0 = 1 + 20 * w;
end

function s = nonlin_S(u, alpha)
s = -u .* (1 - u) + u ./ (alpha + u);
end

%% Solve B2

% Mesh and assemble matrices
geometry = @circleg ;
hmax = 1/40;
[p , e , t] = initmesh(geometry ,'hmax' , hmax);
[A,M,b] = assema(p,t, 1, 1, 0);


% Population count 
Population_Rate = zeros(size(tvec));
% Solve over time with crank-nicolson + fixed-point
u = initial_condition(p(1,:), p(2,:))'; 
for time = tvec(2:end)
    % % fixed-point iteration function handle

    F = @(u_next) ((M/dt + (d1/2)*A)) \ (((M/dt - (d1/2)*A))*u - M * nonlin_S((u_next + u)/2, alpha));
    u_next = fixedPointIteration(F, u, tol, maxIter);

    % Optional: plot intermediate solutions at certain time steps
    if mod(time, 1) < dt
        plotPDE(p,e,t,u, time);
    end

    u = u_next; % Update for next time step
end
%% Plot solution
function plotPDE(p,e,t,u,time)
figure;
pdeplot(p,e,t,'XYData',u,'ZData',u,'Mesh','on');
title('Solution \psi at t = ' + string(time));
xlabel('x');
ylabel('y');
zlabel('\psi');
axis equal;
colorbar;
end
