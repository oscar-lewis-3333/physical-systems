import numpy as np
from scipy.linalg import eigh
from scipy.integrate import simpson
from scipy.interpolate import interp1d

#we start with the solution in 1D

def solve_tise_1d(potential, x, V0, hbar=1, m=1): #we choose hbar=m=1 for scaling reasons (can non-dimensionalise TISE/TDSE). given their real values we would encounter numerical approximation innaccuracies due to how small they are
    N = len(x)
    dx = x[1] - x[0]
    #we construct the hamiltonian matrix for this problem inside (discretise the problem). V = 0 inside 
    main_diag = (hbar**2) / (m * dx**2) * np.ones(N)
    off_diag = -0.5 * hbar**2 / (m * dx**2) * np.ones(N-1)
    H = np.diag(main_diag) + np.diag(off_diag, k=1) + np.diag(off_diag, k=-1)
    V = potential(x)
    H += np.diag(V)


    energies, wf = eigh(H) #in QM, energy states of a particle are the eigenvalues of H, wavefunction at that point in time is the corresponding eigenvectors
    #obviously the eigenvalues are functions here, but they are still vectors in the correctly defined Hilbert space (L^2([a,b]))

    for i in range(len(energies)):
        norm = simpson(wf[:, i]**2, x) #gives the (approximate) norm^2 of the wavefunction with the first i components. 
        wf[:, i] /= np.sqrt(norm)

    #we only keep the bound states

    bound = energies < V0
    return energies[bound], wf[:]

def potential_well_1d(x, V0 = 10.0, a=2.0):
    #V=0 inside, V=V0 outside
    V = np.full_like(x, V0)

    V[np.abs(x) < a] = 0

    return V    


def probability_outside_1d(x_grid, phi_1d, E_1d, a):
    outside_mark = np.abs(x_grid) > a
    dx = x_grid[1] - x_grid[0]
    outside_prob = []
    for i in range(len(E_1d)):
        prob_out = simpson(phi_1d[:, i] **2 * outside_mark, x=x_grid)
        outside_prob.append(prob_out)
    return outside_prob

#we now look at the real part of the time evolved wavefunction
def real_time_evolved_1d(phi_n, E_n, x_grid, t_max, N_t, hbar=1.0):
    #phi_n, E_n, x_grid, t_max, hbar as already defined/logical. N_t number of timesteps. Return meshgrid arrays and real part of time evolved 
    t_arr = np.linspace(0, t_max, N_t)
    #use solution psi_n(x, t) = phi_n(x) exp(-i E_n t/hbar)
    psi_t = phi_n[:, np.newaxis] * np.exp(-1j * E_n * t_arr /hbar)
    Re_psi = np.real(psi_t).T
    X, T = np.meshgrid(x_grid, t_arr)
    return X, T, Re_psi

def superposition_prob(phi_1, E_1, phi_2, E_2, x_grid, t_max=10, N_t=300, hbar=1):

    t_arr = np.linspace(0, t_max, N_t)
    phi_t = (phi_1[:, np.newaxis] * np.exp(-1j * E_1 * t_arr / hbar) + phi_2[:, np.newaxis] * np.exp(-1j * E_2 * t_arr / hbar))
    prob = np.abs(phi_t)**2 #prob density
    return t_arr, prob.T

def solve_3d_well(a, b, c, V0, N_1d, N_3d):
    #we solve by splitting the problem into 3 1D problems and combining back together.

    hx, hy, hz = a/2, b/2, c/2

    pad = 3.0

    x1 = np.linspace(-hx - pad, hx + pad, N_1d)
    y1 = np.linspace(-hy - pad, hy + pad, N_1d)
    z1 = np.linspace(-hz - pad, hz + pad, N_1d)

    Vx = lambda x: potential_well_1d(x, V0, hx)
    Vy = lambda y: potential_well_1d(y, V0, hy)
    Vz = lambda z: potential_well_1d(z, V0, hz)

    E_x, phi_x = solve_tise_1d(Vx, x1, V0)
    E_y, phi_y = solve_tise_1d(Vy, y1, V0)
    E_z, phi_z = solve_tise_1d(Vz, z1, V0)

    phi_x_int = [interp1d(x1, phi_x[:, i], kind='cubic', bounds_error=False, fill_value=0.0) for i in range(len(E_x))]
    phi_y_int = [interp1d(x1, phi_y[:, j], kind='cubic', bounds_error=False, fill_value=0.0) for j in range(len(E_y))]
    phi_z_int = [interp1d(x1, phi_z[:, k], kind='cubic', bounds_error=False, fill_value=0.0) for k in range(len(E_z))]   

    x3 = np.linspace(x1.min(), x1.max(), N_3d)
    y3 = np.linspace(y1.min(), y1.max(), N_3d)
    z3 = np.linspace(z1.min(), z1.max(), N_3d)

    X3, Y3, Z3 = np.meshgrid(x3, y3, z3, indexing='ij')

    return X3, Y3, Z3, (E_x, E_y, E_z), (phi_x_int, phi_y_int, phi_z_int)

def compute_3d_state_seperable(X3, Y3, Z3, interps, indices):
    phi_x_int, phi_y_int, phi_z_int = interps
    nx, ny, nz = indices
    phi_3d = phi_x_int[nx](X3) * phi_y_int[ny](Y3) * phi_z_int[nz](Z3)
    prob_3d = np.abs(phi_3d) **2
    return phi_3d, prob_3d

