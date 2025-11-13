function B=load_vector(x, func, u_a, u_b)
%
% Returns the assembled load vector b.
% Input is a vector x of node coords, and boundary values u_a, u_b.
%
N = length(x) - 1;
B = zeros(N+1, 1);
for i = 1:N
    h = x(i+1) - x(i);
    n = [i i+1];
    B(n) = B(n) + [func(x(i)); func(x(i+1))]*h/2; % Trapezoidal rule
end
B(1) = u_a;
B(N+1) = u_b;
end
