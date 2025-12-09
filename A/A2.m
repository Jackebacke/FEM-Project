clear; clc; close all;
%% Equation params
a = -1; b = 1; % end points of interval
u_a = 0; u_b = 0; % boundary value at x=a and x=b
delta = 0.1;

%% FEM params
TOL = 1e-6;
max_nodes = 1e4;
nodes_start = 12;
lambda = 0.9;

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
%% Main script

% Initialize variables
x = linspace(a, b, nodes_start); % initial mesh coords
NODES = [length(x)]; % History of #nodes
ETA = [1]; % Error estimator history
iteration = 0;

while ETA(end) > TOL && length(x) < max_nodes
    fprintf("Iteration %d: Number of nodes = %d, Error estimator = %.5e\n", iteration, NODES(end), ETA(end));
    N = length(x); % number of elements

    A = stiffness_matrix(x);
    B = load_vector(x, @(x) (f(x)/delta), u_a, u_b);
    u_h = A\B; % solution

    M = mass_matrix(x);
    Zeta = -M \ A*u_h; % Second derivative approximation at nodes

    % Compute error estimator

    eta2 = zeros(N,1);
    for i = 1:N-1
        h = x(i+1) - x(i);
        eta2(i) = (h^3 / 2) * (abs(f(x(i)))^2 + abs(f(x(i+1)))^2); % no laplacian
        % eta2(i) = (h^3 / 2) * (abs(f(x(i)) + delta * Zeta(i))^2 + abs(f(x(i+1)) + delta * Zeta(i+1))^2); % Trapezoidal rule
    end


    if sum(eta2) < TOL
        disp('Desired tolerance reached. Stopping refinement.');
        break;
    end

    % Mesh refinement
    for i = 1: length(eta2)
        if eta2(i) > lambda * max(eta2) && length(x) < max_nodes
            x = [x (x(i+1) + x(i))/2];
        end
    end
    x = sort(x);

    % Handle if max_nodes is reached
    if length(x) >= max_nodes
        display('Maximum number of nodes reached. Stopping refinement.');
        N = length(x); % number of elements

        A = stiffness_matrix(x);
        B = load_vector(x, @(x) (f(x)/delta), u_a, u_b);
        u_h = A\B; % solution

        M = mass_matrix(x);
        Zeta = -M \ A*u_h; % Second derivative approximation at nodes

        % Compute error estimator

        eta2 = zeros(N,1);
        for i = 1:N-1
            h = x(i+1) - x(i);
            eta2(i) = (h^3 / 2) * (abs(f(x(i)) + delta * Zeta(i))^2 + abs(f(x(i+1)) + delta * Zeta(i+1))^2); % Trapezoidal rule
        end
    end

    NODES = [NODES length(x)]; % Update history of #nodes
    ETA = [ETA sum(eta2)]; % Update error estimator history
    iteration = iteration + 1;
end

disp('---------------- Final solution computed. -----------------');
disp(['Total iterations: ' num2str(iteration)]);
disp(['Final number of nodes: ' num2str(NODES(end))]);
disp(['Final error estimator: ' num2str(ETA(end))]);


%% Plots
Residual = f(x)' + delta * Zeta; % Residual at nodes

figure();
plot(x, delta * Zeta, '-o');
hold on;
plot(x, f(x), '-x');
legend('\delta \cdot \zeta (x)', 'f(x)');
xlabel('x');
ylabel('Second derivative value');
title('Second Derivative Approximation');


figure();
% Set figsize for better printing
set(gcf, 'Position', [100, 100, 800, 1000]);
% Subplots
subplot(4,1,1);
plot(x, u_h, '-o');
xlabel('x');
ylabel('u_h(x)');
title('FEM Solution u_h(x)');

subplot(4,1,2);
plot(x, Residual, '-o');
xlabel('x');
ylabel('Residual');
title('Residual at nodes');

subplot(4,1,3);
plot(x, sqrt(eta2), '-o');
xlabel('x');
ylabel('Element residual');
title('Error indicator \eta_i per element');

subplot(4,1,4);
plot(x(2:end), 1./diff(x), '-o');
xlabel('x');
ylabel('Grid size');
title('Grid size distribution');

figure();
loglog(NODES, ETA, 'or-');
hold on;
loglog(NODES, NODES.^-1, 'ob-');

% Compute convergence rate of ETA
logN = log(NODES(2:end));  % skip initial point
logETA = log(ETA(2:end));
p = polyfit(logN, logETA, 1);  % linear fit
rate = p(1);  % convergence rate
const = exp(p(2));
fprintf('Estimated convergence rate of ETA: ETA ~ {%.3f} * N^{%.3f})\n', const, rate);
% Add fitted line to plot
loglog(NODES, const * NODES.^rate, '--k', 'LineWidth', 1.5);

legend('\eta', 'N^{-1} reference', sprintf('Fit: {%.2e} \\cdot N^{%.2f}', const, rate));
xlabel('Number of nodes');
ylabel('Sum of error indicators');
title('Sum of error indicators vs. Number of nodes');


