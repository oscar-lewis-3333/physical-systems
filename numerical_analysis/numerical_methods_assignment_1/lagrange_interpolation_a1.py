import numpy as np

def lagrange (k, x, nodes):
    terms = [(x - nodes[j])/(nodes[k] - nodes[j]) for j in range(len(nodes)) if j!=k]
    return np.prod(terms)

def lagrange_interp(x, nodes, values):
    result = np.zeros_like(x)  # x could be an array of any shape, so we create an array of zeros in the same shape as x in order to add the contributions of each lagrange interpolation to 
    for k in range(len(nodes)):  
        result += values[k] * np.array([lagrange(k, xi, nodes) for xi in x])  #summing together
    return result

def q(x):
    return 1/ (1 + x**2)

def m_value(m):
    nodes = np.linspace(-5, 5, m)
    values = q(nodes)
    return nodes, values
