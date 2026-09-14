import numpy as np

def f(t,y):
    return -(1 + t)* y**2

def exact_solution(t):
    return 2/(t**2 + 2*t + 2)

def evolve(t0, y0, f, Phi, h, N): #my evolve function from the quiz
    u = np.zeros(N+1)
    u[0] = y0
    for n in range(N):
        t_n = t0 + n * h
        u_n = u[n]
        u[n+1] = u_n + h * Phi(t_n, u_n, f, h)
    return u

def last_error(Phi, h, N):
    t0 = 0 #given initial conditions
    y0 = 1
    u = evolve(0, 1, f, Phi, h, N) #approximate solution
    t_N = N*h  #t_N = t0 + N*h = N*h
    y_N = exact_solution(t_N) 
    return np.abs(y_N - u[-1]) #u[-1] is the last element of the array, and hence corresponds to u_N

def EOC(hs, es):
    hs = np.array(hs) #making hs, es into arrays (vectors)
    es = np.array(es)
    eoc = np.zeros(len(hs)) #creating an array for eoc 
    eoc[0] = np.nan #first value doesn't make any sense so we put this suitable value there instead
    for i in range(1, len(hs)):
        eoc[i] = np.log(es[i]/es[i-1])/np.log(hs[i]/hs[i-1]) #by definition of eoc
    return eoc

k_values =np.arange(11)

h_values = [(1.5/(4*(2**k))) for k in k_values]

N_values = [(4 * (2**k)) for k in k_values]

def errors(Phi, h, N):
    return [last_error(Phi, h, N) for h, N in zip(h_values, N_values)]

def Euler(t, y, f, h): #euler increment from the quiz
    return f(t,y)

def Heun(t, y, f, h): #heun increment function from the quiz
    return 0.5 * (f(t,y) + f(t+h, y + h*f(t,y)))
