clear; clc; close all;
%% Geometry
geometry = @circleg ;

%% Functions
a = @(x,y) 1.0;
f = @(x,y) 8.*pi^2 .* sin(2.*pi.*x) .* sin(2.*pi.*y);
u_exact = @(x,y) sin(2.*pi.*x) .* sin(2.*pi.*y);

function EnE = EnergyError(p, A, u_exact, u_h)
% Evaluate exact solution at mesh nodes
disp('fgjjdfgdkfgjdfg');
size(u_h)
u_ex = u_exact(p(1,:), p(2,:))';
size(u_ex)
% Compute error
err = u_ex - u_h;

% Energy norm
EnE = sqrt(err' * A * err);
end

H_MAX = [1/2, 1/4, 1/8, 1/16, 1/32];
errors = [];
for hmax = H_MAX
    [p , e , t] = initmesh(geometry ,'hmax' , hmax);
    I = eye ( length ( p )); % construct the identity matrix
    A = assembleA(p,e,t,a); % stiffness matrix
    b = assembleb(p,e,t,f); % load vector
    A(e(1 ,:),:) = I(e(1 ,:) ,:); % apply Dirichlet BCs strongly
    b(e(1 ,:)) = u_exact(p(1 ,e(1 ,:)) , p(2 ,e(1 ,:))); % apply Dirichlet BCs strongly
    size(b(e(1 ,:)))
    u_h = A \ b;
    errors = [errors, EnergyError(p, A, u_exact, u_h)];
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
