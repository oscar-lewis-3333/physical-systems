import numpy as np

def recursive_f(f, nodes):
    if len(nodes) == 1:
        return f(nodes[0])
    else:
        f_1 = recursive_f(f, nodes[1:]) #nodes[1:] is every node apart from the first one
        f_2 = recursive_f(f, nodes[:-1])
        return (f_1 - f_2)/ (nodes[-1] - nodes[0]) #Using the formula in 2.4

def h(x):
    return np.sin(x)

def qk(f, nodes, k, x):
    if k == 0:
        return f(nodes[0]) #q_0 (x) = f(x_0)
    else:
        q = qk(f, nodes[:-1], k - 1, x) #Defining q_k-1 (x) in order to recursively define q_k (x)
        f_recursive = recursive_f(f, nodes[:k+1]) #Nodes[:k+1] slices off the nodes after x_0,...,x_k giving the correct set
        terms = np.prod([(x - nodes[j]) for j in range(k)]) #(x - x_0)...(x - x_k-1)
        return q + f_recursive * terms

def j(x):
    return np.arctan(x)
