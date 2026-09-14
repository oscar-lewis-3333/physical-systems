import numpy as np

def exact_solution(alpha, t):
    return np.exp(alpha * t) #solving the ODE

def forward_euler(t_s, t_e, N, alpha):
    h = (t_e - t_s)/N #equally split time-steps
    u = np.zeros(N+1) #creating an array of zeros in order to store our values of u_n
    u[0] = 1  #initial condition
    for n in range(N):
        u[n+1] = u[n] + h*alpha*u[n] #using fn = f(t_n, u_n) = alpha*u_n
    return u

def heun(t_s, t_e, N, alpha):
    h = (t_e - t_s) / N #equal spaced time steps
    u = np.zeros(N+1) #as for the forward euler, we create an array of zeros to store our values of u_n+1
    u[0] = 1 #initial condition
    for n in range(N):
        u[n+1] = u[n] + (h/2)*( alpha*u[n] + alpha*(u[n] + h*alpha*u[n])) #using the formula, fn = alpha*u_n and f(t_n+1, u_n + hf_n) = alpha* (u_n + h*(alpha*u_n))
    return u
