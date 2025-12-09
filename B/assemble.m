function [A,R,b,r] = assemble(p,e,t, a, f, g, gamma)
% Assembles the global stiffness matrix A, boundary mass matrix R,
% load vector b, and boundary load vector r for the problem
% -div(a grad u) = f in Omega
% a du/dn = gamma(g - u) on dOmega
% Input:
%   p, e, t - mesh Data
%   a       - function handle for coefficient a(x,y)
%   f       - function handle for right-hand side f(x,y)
%   g       - function handle for boundary data g(x,y)
%   gamma   - function handle for boundary coefficient gamma(x,y)
% Output:
%   A       - global stiffness matrix
%   R       - boundary mass matrix
%   b       - load vector
%   r       - boundary load vector
N = size(p,2);
A = sparse(N,N);
R = sparse(N,N);
b = zeros(N,1);
r = zeros(N,1);
% assemble stiffness matrix A, and load vector b.
for K = 1:size(t,2);
    nodes = t(1:3,K);
    x = p(1,nodes);
    y = p(2,nodes);
    [area,bK,cK] = HatGradients(x,y);
    xc = mean(x); yc = mean(y); % element centroid

    abar = a(xc,yc); % value of a(x,y) at centroid
    AK = abar*(bK*bK' + cK*cK')*area; % element stiffness matrix
    bK = f(xc,yc)*area/3;
    A(nodes,nodes) = A(nodes,nodes) + AK;
    b(nodes) = b(nodes) + bK;
end

% assemble boundary mass matrix R, and the vector r.
for E = 1:size(e,2)
    nodes = e(1:2,E);
    x = p(1,nodes);
    y = p(2,nodes);
    length_E = sqrt((x(1)-x(2))^2+(y(1)-y(2))^2);
    R(nodes,nodes) = R(nodes,nodes) + gamma(x,y)*length_E*[2 1;1 2]/6;
    r(nodes) = r(nodes) + gamma(x,y)*length_E*g(x,y)'/2;
end
end



