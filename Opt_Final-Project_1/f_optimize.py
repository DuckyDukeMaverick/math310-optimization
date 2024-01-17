# Tran Khue Tu & Nguyen Viet Duc
# MATH301: Optimization 

import numpy as np
import time
from time import process_time
NORM = np.linalg.norm

# Weak wolfe line search
def EBLS(f, gradf, current_point, direction):
    c_1 = 1e-6
    c_2 = 0.8
    alpha = 1
    L = 0
    U = 10
    x = np.array(current_point)
    d = np.array(direction)
    gradf_curr = np.array(gradf(x))

    # Terminate when cannot find alpha satisfy 
    for i in range(10000):
        if f(x + alpha*d) > f(x) + c_1*alpha*np.dot(gradf_curr, x):
            U = alpha
            alpha = (U + L)/2
        elif np.dot(np.array(gradf(x + alpha*d)), d) < c_2*alpha*np.dot(gradf_curr, d):
            L = alpha
            if U >= 10:
                alpha = 2*L
            else:
                alpha = (L + U)/2
        else:
            return alpha
    return None


# Conjugate Gradient for NonLinear: Polak Ribiere 
def f_optimize(f,df,x0,time_limit):
    # Parameters 
    n = len(x0)
    start_time = time.time()
    max_iterations = 10000
    grad_tol = 1e-6

    # Initialize storage for n variables, starting with x0 
    xs = ["x%d" % x for x in range(n)]
    for i in range(n):
      xs[i] = [x0[i]]

    # Initialize the descent direction
    D = df(x0)
    delta = -1 * D 

    current_time = time.time()
    while ((current_time - start_time < time_limit*0.98)):
        x = x0
        alpha = EBLS(f=f, gradf=df, current_point = x, direction=delta)

        if alpha!=None:
            X = x0+ alpha*delta # Newly updated experimental point
        if NORM(df(X)) < grad_tol:
            for i in range(n):
              xs[i] += [X[i], ]
            return X # Return the results
        else:
            x0 = X
            d = D # Gradient of the preceding experimental point
            D = df(x0) # Gradient of the current experimental point
            chi = np.array(D-d).dot(D)/NORM(d)**2
            chi = max(0, chi) # Line (16) of the Polak-Ribiere Algorithm
            delta = -D + chi*delta # Newly updated direction
            for i in range(n):
              xs[i] += [X[i], ]

        current_time = time.time()

    return X

