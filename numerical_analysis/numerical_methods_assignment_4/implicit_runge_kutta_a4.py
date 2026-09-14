import numpy as np
from scipy.optimize import root_scalar

def DIRK(t, u, f, h, A, b, c):
    s = len(b)
    K = np.zeros(s)
    for i in range(s):
        def F(x):
            t_i = t + c[i]*h #formula for the first entry of f in F(x) formula
            u_i = u + h* np.sum(A[i][j]*K[j] for j in range(i)) + h* A[i][i]*x
            return x - f(t_i, u_i)
        K[i] = root_scalar(F, x0 = f(t,u), x1 = 2*f(t,u)).root
    return np.sum(b[i] * K[i] for i in range(s))

def implicit_midpoint(t, u, f, h):
    A = [[1/2]]
    b = [1]
    c = [1/2]
    return DIRK(t, u, f, h, A, b, c)

def crouzeix(t, u, f, h):
    A = [[0.5 + (np.sqrt(3)/6), 0], [-np.sqrt(3)/3, 0.5 + (np.sqrt(3)/6)]]
    b = [0.5, 0.5]
    c = [0.5 + np.sqrt(3)/6, 0.5 - np.sqrt(3)/6]
    return DIRK(t, u, f, h, A, b, c)
