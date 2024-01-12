function [p, u, t_elapsed] = interior_point(A, b, c)
time = tic;
Aeq = [];
beq = [];
lb = zeros(size(c));
ub = ones(size(c)) * abs(max(A, [], 'all'));
options = optimoptions('linprog','Algorithm','interior-point', 'MaxIter',10000);

[x, fval] = linprog(c,A,b,Aeq,beq,lb,ub,options);

u = 1/fval;
p = x*u;
t_elapsed = toc(time);