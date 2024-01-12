% newton's method
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

% function part a
function f = hw5_1a_f(x)
f = x^4-1;
function v = hw5_1a_gf(x)
v = 4*x^3;
function v = hw5_1a_inv_ggf(x)
v = 1/(12*x^2);

% function part b
function f = hw5_1b_f(x)
f = (10*x(1)^2+ x(2)^2)/2;
function gf = hw5_1b_gf(x)
gf = transpose([10*x(1) x(2)]);
function inv_ggf = hw5_1b_inv_ggf(x)
inv_ggf = [1/10 0; 0 1];

% function part c
function f = hw5_1c_f(x)
f = log(exp(x) + exp(-x));
function gf = hw5_1c_gf(x)
gf = (exp(x) - exp(-x))/(exp(x) + exp(-x));
function inv_ggf = hw5_1c_inv_ggf(x)
inv_ggf = ((exp(x) + exp(-x))^2)/4;

% function part d
function f = hw5_1d_f(x)
f = -log(x) + x;
function gf = hw5_1d_gf(x)
gf = -1/x + 1;
function inv_ggf = hw5_1d_inv_ggf(x)
inv_ggf = (x^2);

% testing with different functions
% part a
[sol_a, f_path_a, path_a] = new_meth(@ hw5_1a_f, @ hw5_1a_gf,...
    @ hw5_1a_inv_ggf, 4, 1e-50, 50, 500);
% part b
[sol_b1, f_path_b1, path_b1] = new_meth(@ hw5_1b_f, @ hw5_1b_gf,...
    @ hw5_1b_inv_ggf, [0 ; 0], 1e-50, 50, 500);
[sol_b2, f_path_b2, path_b2] = new_meth(@ hw5_1b_f, @ hw5_1b_gf,...
    @ hw5_1b_inv_ggf, [10 ; 10], 1e-50, 50, 500);
% part c
[sol_c1, f_path_c1, path_c1] = new_meth(@ hw5_1c_f, @ hw5_1c_gf,...
    @ hw5_1c_inv_ggf, 1, 1e-50, 50, 500);
[sol_c2, f_path_c2, path_c2] = new_meth(@ hw5_1c_f, @ hw5_1c_gf,...
    @ hw5_1c_inv_ggf, 1.1, 1e-50, 50, 500);

% part d
[sol_d, f_path_d, path_d] = new_meth(@ hw5_1d_f, @ hw5_1d_gf,...
    @ hw5_1d_inv_ggf, 3, 1e-50, 50, 500);