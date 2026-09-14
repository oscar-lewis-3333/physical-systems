import numpy as np

def rotating_hoop(t, y, R, g, Omega):
    theta, theta_dot = y

    theta_ddot = np.sin(theta) * (Omega **2 * np.cos(theta) - g/R)

    return [theta_dot, theta_ddot]

def solution(R, g, Omega, y0, t, rtol=1e-9, atol=1e-11):
    from scipy.integrate import solve_ivp
    sol = solve_ivp(rotating_hoop, t_span=(t[0], t[-1]), y0=y0, t_eval=t, args=(R, g, Omega), method="DOP853", rtol=rtol, atol=atol)

    if not sol.success:
        raise RuntimeError(f"Numerical integration failed: {sol.message}")

    return sol

def energy(y, m, R, g, Omega):
    theta, theta_dot = y[0], y[1]
    E = 0.5 * m * R**2 * theta_dot**2 - 0.5 * m * R**2 * Omega**2 * (np.sin(theta))**2 - m*g*R*np.cos(theta) 
    return E

def dimensionless_potential(theta, Omega, Omega_c):
    ratio = Omega / Omega_c

    return (1 - np.cos(theta) - 0.5 * ratio**2 * np.sin(theta)**2)

def dimensionless_inital_energy(y, Omega, Omega_c):
    theta_0 = y[0]
    theta_dot_0 = y[1]

    energy_0 = theta_dot_0**2 / (2 * Omega_c**2) + dimensionless_potential(theta_0, Omega, Omega_c)
    return energy_0