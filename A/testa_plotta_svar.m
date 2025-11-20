clear; clc; close all;
%% Equation parameters
a = -1; b = 1; % end points of interval
u_a = 0; u_b = 0; % boundary value at x=a and x=b
delta = 0.1;

%% FEM parameters

TOL = 1e-3;
max_nodes = 1e4;
N_start = 10000;
lambda = 0.9;
h = (b-a)/N_start; % initial mesh size
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


A = stiffness_matrix(x);
B = load_vector(x, @(x) (f(x)/delta), u_a, u_b);
M = mass_matrix(x);

u_h = A\B; % solution
Zeta = -M \ A*u_h;

figure();
plot(x, u_h, '-o');
hold on;
plot(x, Zeta, '-x');
plot(x, f(x), '-s');
xlabel('x');
ylabel('u_h(x), Zeta(x), f(x)');
title('FEM Solution u_h(x), Zeta(x)');
legend('u_h(x)', 'Zeta(x)', 'f(x)');




