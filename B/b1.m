clear; clc; close all;
geometry = @circleg ;
hmax = 1/5;
[p , e , t] = initmesh(geometry ,'hmax' , hmax);

kappa = @(x,y) 1.0;
f = @(x,y) 8.*pi^2 .* sin(2.*pi.*x) .* sin(2.*pi.*y);
gamma = @(x,y) 1; % Large gamma to approximate Dirichlet BCs
g = @(x,y) (2.*x.*pi.*cos(2.*pi.*x).*sin(2.*pi.*y) + 2.*y.*pi.*sin(2.*pi.*x).*cos(2.*pi.*y)) ./ (sqrt(x.^2 + y.^2)); % Neumann BC on boundary

[A , R , b , r ] = assemble( p , e , t , kappa , f , g , gamma)

u_h = (A) \ (b + r);

%% Plotting the solution
figure;
pdeplot(p,e,t,'XYData',u_h,'ZData',u_h,'Mesh','on');
title('FEM Solution: -\Delta u = f with Neumann BCs');
xlabel('x');
ylabel('y');
zlabel('u_h');
axis equal;
colorbar;