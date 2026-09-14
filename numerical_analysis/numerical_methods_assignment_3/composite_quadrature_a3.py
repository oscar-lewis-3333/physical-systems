import numpy as np

def scientific (a, precision =3): #helper code given in question
    expo =int(np.log(abs(a))/np.log (10))
    expo = expo if expo >0 else expo -1
    manti = np.round(float(a)/10** expo, precision)
    if expo < -1 or expo > 1:
        return '{}E{}'.format(manti, expo)
    else:
        return '{}'.format(np.round (a, precision))

def composite_trapezoidal(f, a, b, n):
    h = (b - a) / n #equally split intervals
    total = 0 #starting a sum
    for j in range(1, n + 1):
        x_i = a + (j - 1)* h #x_j-1
        x_j = a + j*h
        total += 0.5 * (x_j - x_i) * (f(x_i) + f(x_j)) #using the trapezoidal formula for each subinterval and summing
    return total

def composite_simpson(f, a, b, n):
    h = (b - a) / n  #equally split intervals
    total = 0 #starting a sum
    for j in range(1, n + 1):
        x_i = a + (j - 1)* h #x_j-1
        x_j = a + j*h
        total += ((x_j - x_i) / 6) * (f(x_i) + 4 * f(0.5*(x_i + x_j)) + f(x_j)) #using simpson formula for each subinterval and summing
    return total

def f(x):
    return np.tan(x)

def real_value(f, a, b):
    return -np.log(np.cos(b)) + np.log(np.cos(a))  #integrating tanx on [a,b] gives the formula
