import numpy as np

def x_plus(p,q):
    result = -1*p + np.sqrt(p**2 + q)
    return result

def x_minus(p,q):
    result = -1*p - np.sqrt(p**2 + q)
    return result

def R_plus(p, dp, q):
    result = ((abs(x_plus(p + dp,q) - x_plus(p,q)))/(abs(x_plus(p,q))))/(abs(dp)/abs(p))
    return result

def R_minus(p, dp, q):
    result = (abs(x_minus(p + dp,q) - x_minus(p,q)))/abs(x_minus(p,q))/(abs(dp)/abs(p))
    return result

def K(p):
    num = abs(p)
    den = np.sqrt(p**2 + 1)
    result = num/den
    return result

def R_plus2(p, dq, q):
    result = ((abs(x_plus(p,q + dq) - x_plus(p,q)))/(abs(x_plus(p,q))))/(abs(dq)/abs(q))
    return result

def R_minus2(p, dq, q):
    result = (abs(x_minus(p,q + dq) - x_minus(p,q)))/abs(x_minus(p,q))/(abs(dq)/abs(q))
    return result

def K_plus(q):
    num = abs(1 + np.sqrt(q+1))
    den = 2*np.sqrt(q+1)
    result = num/den
    return result

def K_minus(q):
    num = abs(1 - np.sqrt(q + 1))
    den = 2*np.sqrt(q+1)
    result = num/den
    return result
