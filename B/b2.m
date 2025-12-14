clear; clc; close all;
%% Params
T = 2;
dt = 0.01;
d1 = 0.1;
alpha = 10;

geometry = @circleg ;

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

function pr = population(u, p, t)
pr = 0;
for K = 1:size(t,2)
    nodes = t(1:3,K);
    x = p(1,nodes);
    y = p(2,nodes);
    area = polyarea(x,y);
    u_avg = mean(u(nodes));
    pr = pr + u_avg * area;
end
end

%% Solve B2
Hs = [1/5, 1/20, 1/40];
pop_timelines = cell(length(Hs),1);
for hmax = Hs
% Mesh and assemble matrices
[p , e , t] = initmesh(geometry ,'hmax' , hmax);
[A,M,b] = assema(p,t, 1, 1, 0);

% Population count
Pop = zeros(size(tvec));

% Solve over time with crank-nicolson + fixed-point
u = initial_condition(p(1,:), p(2,:))';
Pop(1) = population(u, p, t);
for time = tvec
    % % fixed-point iteration function handle
    F = @(u_next) ((M/dt + (d1/2)*A)) \ (((M/dt - (d1/2)*A))*u - M * nonlin_S((u_next + u)/2, alpha));
    u_next = fixedPointIteration(F, u, tol, maxIter);

    % Optional: plot intermediate solutions at certain time steps
    if time == 0 || time == T
        plotPDE(p,e,t,u, time, hmax);
    end

    u = u_next; % Update for next time step
    Pop(find(tvec == time)) = population(u, p, t);
end
pop_timelines{find(Hs == hmax)} = Pop;
end

%% Plot solution
function plotPDE(p,e,t,u,time, hmax)
figure;
pdeplot(p,e,t,'XYData',u,'ZData',u,'Mesh','on');
title('Solution v(x,t) at T = ' + string(time) + ', hmax = ' + string(hmax));
xlabel('x');
ylabel('y');
zlabel('v(x,t)');
colorbar;
name = sprintf('B2_solution_t_%.2f_hmax_%.4f.pdf', time, hmax);
folder = './Project/B/Figures';
savefig_report(name, folder, 'Width', 7, 'Height', 9);
end

%% Plot population rate over time
figure;
for i = 1:length(Hs)
    hold on;
    Pop = pop_timelines{i};
    plot(tvec, Pop, '-', 'DisplayName', 'hmax = ' + string(Hs(i)));
end
title('Population over time');
xlabel('Time');
ylabel('Population');
grid on;
legend('show');
savefig_report('B2_population_over_time.pdf', './Project/B/Figures', 'Width', 7, 'Height', 9);
