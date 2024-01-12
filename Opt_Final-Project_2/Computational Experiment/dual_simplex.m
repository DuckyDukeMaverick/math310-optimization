function [p, u, t_elapsed] = dual_simplex(A, b, c)
time = tic;% star measuring time 
lb = zeros(size(c));
ub = ones(size(c)) * abs(max(A, [], 'all'));

[x, fval] = linprog(c,A,b, [], [], lb, ub);
u = 1/fval;
p = x*u;
t_elapsed = toc(time);