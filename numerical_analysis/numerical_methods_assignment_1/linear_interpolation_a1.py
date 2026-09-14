import numpy as np

def linear_interpolate(f, x_0, x_1, x):
    result = f(x_0) + ((f(x_1) - f(x_0))*(x - x_0))/(x_1 - x_0)
    return result

def g(x):
    return 2* np.sin(2*x)

def p_1(x):
    result = linear_interpolate(g, 0, 1, x)
    return result

def h(x):
    return np.sin(x)

def p_2(x):
    return linear_interpolate(h, 0, 1, x)

def interpolation_error(f, x_0, x_1, x):
    return f(x) - linear_interpolate(f, x_0, x_1, x)

def j(x):
    return np.log(x)

def r(x):
    return 1/(1+x**2)

def max_error(f, x_0, x_1):
    return max(abs(interpolation_error(f, x_0, x_1, x)) for x in np.linspace(x_0, x_1, 101))
