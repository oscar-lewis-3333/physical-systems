import numpy as np

def explicitRK(t, u, f, h, A, b, c):
    s = len(b) #by definition of the Butcher Tableau
    K = np.zeros(s) #creating an array to store values
    for i in range(s):
        t_i = t + c[i]*h #first function entry in definition of K_i
        u_i = u + h * np.sum(A[i][j] * K[j] for j in range(i+1)) #second function entry in definition of K_i, but we only sum to i since the matrix is lower-triangular, hence every other entry is zero.
        K[i] = f(t_i, u_i)
    Phi = np.sum(K[i] * b[i] for i in range(s))
    return Phi

def forwardeuler(t, u, f, h):
    A = [[0]] 
    b = [1]
    c = [0]
    return explicitRK(t, u, f, h, A, b, c)

def explicit_midpoint(t, u, f, h):
    A = [[0, 0], [1/2, 0]]
    b = [0, 1]
    c = [0, 1/2]
    return explicitRK(t, u, f, h, A, b, c)

def RK3 (t, u, f, h):
    A = [[0, 0, 0], [1/2, 0, 0], [-1, 2, 0]]
    b = [1/6, 2/3, 1/6]
    c = [0, 1/2, 1]
    return explicitRK(t, u, f, h, A, b, c)
