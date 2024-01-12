for i = 1:1:1000
    m = 100;
    n = 100;
     
    A = randi(100, m,n) * -1;
    b = ones(1, m) * -1;
    c = ones(1, n);
    
    [p_dsimplex, u_dsimplex, t_elapsed_dsimplex] = dual_simplex(A, b, c);
    [p_interior, u_interior, t_elapsed_interior] = interior_point(A, b, c);

    pts(i, 1) = t_elapsed_dsimplex;
    pts(i, 2) = t_elapsed_interior;
    pts(i, 3) = abs(u_dsimplex - u_interior);
end 

boxplot([pts(:, 1), pts(:, 2)], 'Notch','on','Labels',{'Simplex','Interior Point'})
title({'Compare Runtime on 100x100 payoff matrix', 'using Dual Simplex and Interior Point'})
ylabel('Time (seconds)')
fontname('Times New Roman')