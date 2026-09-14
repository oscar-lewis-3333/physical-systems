import numpy as np

def G(f, a, b):
    term1 = (5/9)*f(0.5*(a + b) + (np.sqrt(3/5) / 2)*(b-a))
    term2 = (8/9) * f(0.5*(b+a))
    term3 = (5/9)*f(0.5*(a + b) - (np.sqrt(3/5) / 2)*(b-a))
    return 0.5*(b-a) * (term1 + term2 + term3)

def F(x):
    return np.cos(x)

def exact_integral(x): 
    return np.sin(x) #By integrating and noticing we are using the interval [0,h]
