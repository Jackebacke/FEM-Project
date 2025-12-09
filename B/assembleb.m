function b = assemble(p,e,t,f)

N = size(p,2);
b = zeros(N,1);
% assemble load vector b.
for K = 1:size(t,2);
    nodes = t(1:3,K);
    x = p(1,nodes);
    y = p(2,nodes);
    [area,bK,cK] = HatGradients(x,y);
    xc = mean(x); yc = mean(y); % element centroid
    bK = f(xc,yc)*area/3;
    b(nodes) = b(nodes) + bK;
end
end