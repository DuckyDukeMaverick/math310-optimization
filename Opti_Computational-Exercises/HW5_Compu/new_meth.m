function [sol, f_path, path] = new_meth(f, gf, inv_ggf, x_0, g_tol, max_time, max_iter)

% initializing the model

x_size = size(x_0);
dim = max(x_size(1), x_size(2));

if x_size(1) > x_size(2)    %(standardizing -> storing in x_path) starting point
x = x_0;
else
x = transpose(x_0);
end

time = tic;
x_path = zeros(dim, max_iter); %initialize convergence path
f_val = zeros(1,max_iter); % initialize solution values
iter = 0; %initialize iteration
t_elapsed = toc(time);

%Newton Method

% (norm(gf(x)) >= g_tol) - I remove this condition so that the algorithms
% can be run in 500 iterations
while (t_elapsed <= max_time) && (iter < max_iter)
    iter = iter + 1;
    for i = 1:dim
        x_path(i,iter) = x(i);
    end
    f_val(iter) = f(x);
    p = inv_ggf(x)*gf(x);
    x = x - p;
    t_elapsed = toc(time);
end

%Output results
sol = x;
f_path = f_val(1:iter);
path = x_path(:, 1:iter);
