[X,Y] = meshgrid(-50:0.1:50, -50:0.1:50);
Z = ((X.^2).*10 + (Y.^2))/2;
f = @(x) (10*x(1).^2 + x(2).^2)/2;

gradf =@(x) [10*x(1);x(2)];                

% Choose an initial guess or starting model point, iterations, and step size.; you can choose any guess and play with the initial guess, step size and iterations.
intialGuess = [0.5;20];
iterations = 10;
stepSize = 5;

recordGuesses = [intialGuess];
nextGuess = intialGuess;
for i = 1: iterations
    % we are going in a negative direction (Gradient Descent), so there will be -ve from the previous guess, not +ve
    nextGuess = nextGuess - stepSize*gradf(nextGuess);
    % record the guesses
    recordGuesses = [recordGuesses,nextGuess];

end
contour(X,Y,Z,50)
hold on;
plot(recordGuesses(1,:),recordGuesses(2,:),'o');
hold off;