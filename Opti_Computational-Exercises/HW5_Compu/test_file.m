% testing with different functions
clc; clearvars
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