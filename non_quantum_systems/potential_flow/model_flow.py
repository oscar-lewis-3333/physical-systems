import numpy as np

def uniform_complex_potential(z, U, beta=0.0):
    if U <= 0:
        raise ValueError("Speed must be positive")

    return U * np.exp(-1j * beta)*z

def source_complex_potential(z, Q):
    #issues when plotting before. np.ma.log applied real-valued domain mask before looking at complex, so incorrectly removes Re(z) <0
    z_values = np.asarray(np.ma.getdata(z), dtype=complex)
    with np.errstate(divide="ignore", invalid="ignore"):
        potential = (Q / (2 * np.pi)) * np.log(z_values)

    if np.ma.isMaskedArray(z):
        return np.ma.array(potential, mask=np.ma.getmaskarray(z))
    return potential

def source_complex_velocity(z, Q):
    return Q/(2*np.pi *z) #diff complex potential

def vortex_complex_potential(z, Gamma):
    z_values = np.asarray(np.ma.getdata(z), dtype=complex)
    with np.errstate(divide="ignore", invalid="ignore"):
        potential = (-Gamma * 1j / (2 * np.pi)) * np.log(z_values)
    
    if np.ma.isMaskedArray(z):
        return np.ma.array(potential, mask=np.ma.getmaskarray(z))
    return potential

def vortex_complex_velocity(z, Gamma):
    return -1j*Gamma/(2 * np.pi * z) #differentiating potential wrt to z. just to get arrows on the plot for direction

def dipole_complex_potential(z, D):
    z_values = np.asarray(np.ma.getdata(z), dtype=complex)
    with np.errstate(divide="ignore", invalid="ignore"):
        potential = D / (2 *np.pi * z_values)
        
    if np.ma.isMaskedArray(z):
        return np.ma.array(potential, mask=np.ma.getmaskarray(z))
    return potential

def dipole_complex_velocity(z, D):
    return - D / (2 * np.pi * z**2)
