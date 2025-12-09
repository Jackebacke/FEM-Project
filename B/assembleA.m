function A = assembleA(p,e,t, a)
N = size(p,2);
A = sparse(N,N);


for K = 1:size(t,2);
    nodes = t(1:3,K);
    x = p(1,nodes);
    y = p(2,nodes);
    [area,bK,cK] = HatGradients(x,y);
    xc = mean(x); yc = mean(y); % element centroid

    abar = a(xc,yc); % value of a(x,y) at centroid
    AK = abar*(bK*bK' + cK*cK')*area; % element stiffness matrix
    A(nodes,nodes) = A(nodes,nodes) + AK;
end
end