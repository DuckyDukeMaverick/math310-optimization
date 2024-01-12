%initialize the model
f = @(x) (10*x(1)^2 + x(2).^2)/2;
grad_f = @(y) [10*y(1) ; y(2) ];
descent_direction = @(z) [-9*z(1); -4*z(2)];
initial = [0.5; 20];
k = 2;
x = initial;

% steepest descent
path_steep = horzcat([0;0], initial);
path_f_steep = horzcat(0,f(initial));
step_size = 1/(sqrt(101));
while abs(path_f_steep(k) - path_f_steep(k-1)) > 0.01
 next_x = x - step_size*grad_f(x);
 path_steep = horzcat(path_steep, next_x);
 path_f_steep = horzcat(path_f_steep, f(next_x));
 x = next_x;
 k = k+1;
end

path_steep(:,1) = [];
path_f_steep(:,1) = [];
iter_steep = k - 1;


%backtracking general gradient descent
k = 2;
path_back = horzcat([0;0], initial);
path_f_back = horzcat(0,f(initial));

while abs(path_f_back(k) - path_f_back(k-1)) > 0.01
 %finding step size
 alpha = 0.5 + randi(1,1)/2;
 alpha_reduction = 0.6 + randi(1,1)/3;
 c = 0.2 + randi(1,1)/2;
 while f(x + alpha*descent_direction(x)) > f(x) + c*alpha*grad_f(x)'*descent_direction(x)
  alpha = alpha*alpha_reduction;
 end
 next_x = x + alpha*descent_direction(x);
 path_back = horzcat(path_back, next_x);
 path_f_back = horzcat(path_f_back, f(next_x));
 x = next_x;
 k = k+1;
end

path_back(:,1) = [];
path_f_back(:,1) = [];
iter_back = k - 1;

%creating grid
[X, Y] = meshgrid(-50:0.05:50, -50:0.05:50);
Z = (10*X.^2 + Y.^2)/2;
  %visualization (steepest descent)
contour(X,Y,Z, linspace(-500,500,100));
hold on
plot(path_steep(1,:), path_steep(2,:), '-o', 'MarkerSize', 5, 'LineWidth', 1);
iter_1 = annotation('textbox',[.15 .25 .2 .2],'String',strcat('Iteration: ', num2str(iter_steep)),'FitBoxToText','on');
obj_val = annotation('textbox',[.15 .2 .2 .2],'String',strcat('Result: ', num2str(path_f_steep(iter_steep))),'FitBoxToText','on');
a = sprintf('Solution: [%.2f, %.2f.]^T', path_steep(:,iter_steep));
solu = annotation('textbox',[.15 .15 .2 .2],'String',a ,'FitBoxToText','on');
hold off
axis([-30, 30, -30, 30]);

    %visualization (backtracking)
contour(X,Y,Z, linspace(-500,500,100));
hold on
plot(path_back(1,:), path_back(2,:), '-o', 'MarkerSize', 5, 'LineWidth', 1);
iter_1 = annotation('textbox',[.15 .25 .2 .2],'String',strcat('Iteration: ', num2str(iter_back)),'FitBoxToText','on');
obj_val = annotation('textbox',[.15 .2 .2 .2],'String',strcat('Result: ', num2str(path_f_back(iter_back))),'FitBoxToText','on');
a = sprintf('Solution: [%.2f, %.2f.]^T', path_back(:,iter_back));
solu = annotation('textbox',[.15 .15 .2 .2],'String',a ,'FitBoxToText','on');
hold off
axis([-30, 30, -30, 30]);

%comparing convergence rate
    plot(1:iter_steep, path_f_steep);
hold on
plot(1:iter_back, path_f_back);
hold off
legend('Steepest Descent', 'Backtracking')