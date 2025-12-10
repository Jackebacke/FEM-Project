clear; clc; close all;
%% Geometry
geometry = @circleg ;

%% Functions
kappa = @(x,y) 1.0;
f = @(x,y) 8.*pi^2 .* sin(2.*pi.*x) .* sin(2.*pi.*y);
gamma = @(x,y) 1; % Large gamma to approximate Dirichlet BCs
g = @(x,y) (2.*x.*pi.*cos(2.*pi.*x).*sin(2.*pi.*y) + 2.*y.*pi.*sin(2.*pi.*x).*cos(2.*pi.*y)) ./ (sqrt(x.^2 + y.^2)); % boundary term

u_exact = @(x,y) sin(2.*pi.*x) .* sin(2.*pi.*y);

H_MAX = [1/2, 1/4, 1/8, 1/16, 1/32];
errors = [];
for hmax = H_MAX
    % Mesh
    [p , e , t] = initmesh(geometry ,'hmax' , hmax);

    % Assemble/Solve system
    [A , R , b , r ] = assemble( p , e , t , kappa , f , g , gamma);
    u_h = (A) \ (b + r);

    % Error computation
    err = u_exact(p(1,:), p(2,:))' - u_h;
    EnE = sqrt(err' * A * err);
    errors = [errors, EnE];
    fprintf('hmax: %.5f, Error: %.5f\n', hmax, errors(end));
end


%% Plotting the solution
figure;
pdeplot(p,e,t,'XYData',u_h,'ZData',u_h,'Mesh','on');
title('FEM Solution: -\Delta u = f with Neumann BCs');
xlabel('x');
ylabel('y');
zlabel('u_h');
axis equal;
colorbar;

figure;
loglog(H_MAX, errors, '-o');
hold on;
xlabel('h_{max}');
ylabel('Energy norm of error');
title('Error vs. Mesh Size');
grid on;
