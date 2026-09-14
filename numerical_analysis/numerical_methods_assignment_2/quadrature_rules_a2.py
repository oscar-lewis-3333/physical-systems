import numpy as np

def midpoint(f, a, b):
    return (b - a) * f(0.5 * (a + b))

def trapezoidal(f, a, b):
    return 0.5*(b - a)*(f(a) + f(b))

def simpson(f, a, b):
    return ((b - a) / 6) *(f(a) + 4*f(0.5*(a + b)) + f(b))

def r(x):
    return np.exp(2*x)

def p(x):
    return 1 / (1 + x)

def exact_value(h): #From integrating the function
    return np.log(1 + h)
