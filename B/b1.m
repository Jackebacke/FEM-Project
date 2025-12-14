clear; clc; close all;

%% Geometry
geometry = @circleg ;

%% Functions
a = @(x,y) 1.0;
f = @(x,y) 8.*pi^2 .* sin(2.*pi.*x) .* sin(2.*pi.*y);
u_exact = @(x,y) sin(2.*pi.*x) .* sin(2.*pi.*y);



H_MAX = 1./(2.^(1:1:5));
errors = zeros(size(H_MAX));
for hmax = H_MAX
    % Mesh
    [p , e , t] = initmesh(geometry ,'hmax' , hmax);

    % Assemble/solve system
    I = eye(length(p)); % construct the identity matrix
    A = assembleA(p,e,t,a); % stiffness matrix
    b = assembleb(p,e,t,f); % load vector
    A(e(1 ,:),:) = I(e(1 ,:) ,:); % apply Dirichlet BCs strongly
    b(e(1 ,:)) = u_exact(p(1 ,e(1 ,:)) , p(2 ,e(1 ,:))); % apply Dirichlet BCs strongly
    u_h = A \ b;

    % Error computation
    err = u_exact(p(1,:), p(2,:))' - u_h;
    EnE = sqrt(err' * A * err);
    errors(find(H_MAX == hmax)) = EnE;

    % Plot first and last solutions
    if hmax == H_MAX(1) || hmax == H_MAX(end)
        plotPDE(p,e,t,u_h,hmax);
    end
end


%% Calculate convergence rates 
% Print results
for i = 1:length(H_MAX)
    fprintf('hmax: %.5f, Error: %.5f\n', H_MAX(i), errors(i));
end

rates = log(errors(1:end-1)./errors(2:end)) ./ log(H_MAX(1:end-1)./H_MAX(2:end));
for i = 1:length(rates)
    fprintf('Convergence rate between hmax %.5f and %.5f: %.2f\n', H_MAX(i), H_MAX(i+1), rates(i));
end

% Linear regression
logH = log(H_MAX(:));
logE = log(errors(:));
coeffs = polyfit(logH, logE, 1);
convergence_rate = coeffs(1);
fprintf('Estimated convergence rate p: %.2f\n', convergence_rate);


%% Plots
function plotPDE(p,e,t,u_h,hmax)
    figure;
    pdeplot(p,e,t,'XYData',u_h,'ZData',u_h,'Mesh','on');
    title('FEM Solution with hmax = ' + string(hmax));
    xlabel('x');
    ylabel('y');
    zlabel('u_h');
    axis equal;
    colorbar;
    savefig_report(sprintf('solution_hmax_%.4f.pdf', hmax), ...
    './Project/B/Figures', 'Width', 7, 'Height', 9);
end


figure;
loglog(H_MAX, errors, '-o', 'DisplayName', 'Computed Errors');
hold on;
loglog(H_MAX, exp(coeffs(2)) * H_MAX.^convergence_rate, '--', 'DisplayName', sprintf('Fit: O(h^{%.2f})', convergence_rate));
loglog(H_MAX, exp(coeffs(2)/2)*H_MAX, '--', 'DisplayName', 'O(h)');
xlabel('h_{max}');
ylabel('Energy norm of error');
legend('show', 'Location', 'best');
title('Error vs. Mesh Size');
grid on;
savefig_report('convergence_2.pdf', './Project/B/Figures', 'Width', 7, 'Height', 9);
