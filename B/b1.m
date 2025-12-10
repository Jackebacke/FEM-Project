clear; clc; close all;

%% Geometry
geometry = @circleg ;

%% Functions
a = @(x,y) 1.0;
f = @(x,y) 8.*pi^2 .* sin(2.*pi.*x) .* sin(2.*pi.*y);
u_exact = @(x,y) sin(2.*pi.*x) .* sin(2.*pi.*y);



H_MAX = [1/2, 1/4, 1/8, 1/16, 1/32];
errors = zeros(size(H_MAX));
for hmax = H_MAX
    % Mesh
    [p , e , t] = initmesh(geometry ,'hmax' , hmax);

    % Assemble/solve system
    I = eye ( length ( p )); % construct the identity matrix
    A = assembleA(p,e,t,a); % stiffness matrix
    b = assembleb(p,e,t,f); % load vector
    A(e(1 ,:),:) = I(e(1 ,:) ,:); % apply Dirichlet BCs strongly
    b(e(1 ,:)) = u_exact(p(1 ,e(1 ,:)) , p(2 ,e(1 ,:))); % apply Dirichlet BCs strongly
    u_h = A \ b;

    % Error computation
    err = u_exact(p(1,:), p(2,:))' - u_h;
    EnE = sqrt(err' * A * err);
    errors(find(H_MAX == hmax)) = EnE;
end

% Print results
for i = 1:length(H_MAX)
    fprintf('hmax: %.5f, Error: %.5f\n', H_MAX(i), errors(i));
end

%% Plots
figure;
pdeplot(p,e,t,'XYData',u_h,'ZData',u_h,'Mesh','on');
title('FEM Solution with hmax = ' + string(hmax));
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
