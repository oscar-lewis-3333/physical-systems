from solving_tunneling import real_time_evolved_1d, potential_well_1d, probability_outside_1d, superposition_prob, solve_3d_well, compute_3d_state_seperable
import ipywidgets as widgets
from IPython.display import display
import numpy as np
import plotly.graph_objects as go
from scipy.integrate import simpson

def plot_for_1d_stationary(a, V0, phi_1d, E_1d, x_grid):
    V_1d = potential_well_1d(x_grid, V0=V0, a=a)

    fig = go.Figure()
    fig.add_trace(go.Scatter3d(x=x_grid, y=V_1d, z=np.zeros_like(x_grid), mode='lines', line=dict(color='black', width=5), name='V(x)'))

    x_boundary = [x_grid.min(), x_grid.max()]
    fig.add_trace(go.Scatter3d( x=x_boundary, y=[V0, V0], z=[0, 0], mode='lines', line=dict(color='red', dash='dash', width=2), name=f'E = V₀ = {V0}'))

    for t in [-a, a]:
        fig.add_trace(go.Scatter3d(x=[t, t], y=[0, V0+1], z=[0, 0], mode='lines', line=dict(color='blue', dash='dash', width=2), name=f'x = {t}'))

    outside_probs = probability_outside_1d(x_grid=x_grid, phi_1d=phi_1d, E_1d=E_1d, a=a)

    for i, E in enumerate(E_1d):
        y_vals = np.full_like(x_grid, E)
        z_vals = 0.6 * phi_1d[:, i]
        fig.add_trace(go.Scatter3d(x=x_grid, y=y_vals, z=z_vals, mode='lines', name=f'n={i}, E={E:.2f}, P_out={outside_probs[i]:.2%}', line=dict(width=3)))


    z_range = 0.6 * np.max(np.abs(phi_1d))
    for v, colour in [(-a, 'blue'),(a, 'blue')]:
        yy = np.linspace(0, V0, 100)
        zz = np.linspace(-z_range, z_range, 100)
        Y, Z = np.meshgrid(yy, zz)
        X = np.full_like(Y, v)
        fig.add_trace(go.Surface(x=X, y=Y, z=Z, colorscale=[[0, f'rgba(0,0,255,0.35)'], [1, f'rgba(0,0,255,0.35)']], showscale=False, name=f'Well edge x={v}'))

    for side, x_range in [('left', [-5, -a]), ('right', [a, 5])]:
        xx = np.linspace(x_range[0], x_range[1], 100)
        yy = np.linspace(0, V0, 100)
        X_rect, Y_rect = np.meshgrid(xx, yy)
        Z_rect = np.full_like(X_rect, 0.0)
        fig.add_trace(go.Surface(x=X_rect, y=Y_rect, z=Z_rect, colorscale=[[0, f'rgba(0,0,255,0.35)'], [1, f'rgba(0,0,255,0.35)']], showscale=False, name='Forbidden Region'))

    fig.update_layout(title='1D bound states with tunnelling tails', scene=dict(xaxis_title='x', yaxis_title='Energy', zaxis_title='φ (scaled)'),legend=dict(x=0.8, y=0.9), width=900, height=600)

    fig.show()


def plot_probability_outside_vs_energy(phi_1d, E_1d, x_grid, a, V0):
#plot probability of the particle being outside the well as a function of energy
    outside_probabilities = probability_outside_1d(x_grid=x_grid, phi_1d=phi_1d, E_1d=E_1d, a=a)
    
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=E_1d, y=outside_probabilities, mode='lines+markers', marker=dict(size=8), name='Numerical bound states', hovertemplate='E = %{x:.3f}<br>P(outside) = %{y:.3%}<extra></extra>'))
    fig.add_vline(x=V0, line_dash='dash', line_color='red', annotation_text='barrier height V₀', annotation_position='top left') #line to indicate potential
    fig.update_layout(title='Probability outside the well versus bound-state energy', xaxis_title='Energy E', yaxis_title='P(|x| > a)', yaxis_tickformat='.0%', width=850, height=500)
    fig.show()

def plot_exponential_tail_comparison(phi_n, E_n, x_grid, a, V0, hbar=1.0, m=1.0, side='right'):
    #compare numerical evanescent tail with claimed probability estimate in notebook
    if E_n >= V0:
        raise ValueError('Exponential decay applies only to a bound state with E < V₀.')
    if side not in {'left', 'right'}:
        raise ValueError("side must be either 'left' or 'right'.")

    if side == 'right':
        mask = x_grid > a
        distance_from_edge = x_grid[mask] - a
    else:
        mask = x_grid < -a
        distance_from_edge = -a - x_grid[mask]

    distances = distance_from_edge #distances from edge
    numerical_amplitude = np.abs(phi_n[mask]) #amplitudes at those distances from edge 
     
    valid = numerical_amplitude > np.finfo(float).tiny #cut out zero to make logarithm well defined, then only chose those values
    distances = distances[valid]
    numerical_amplitude = numerical_amplitude[valid]
    if len(distances) == 0:
        raise ValueError('The supplied grid does not contain points outside the selected well edge.')

    kappa = np.sqrt(2 * m * (V0 - E_n)) / hbar #constant defined in notebook. Classically labelled as kappa
    #start theoretical curve at the first exterior grid point. This isolates the predicted decay rate from any discretisation at the edge.

    reference_amplitude = numerical_amplitude[0] * np.exp(-kappa * (distances - distances[0]))

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=distances, y=numerical_amplitude, mode='lines', name='numerical |φ|'))
    fig.add_trace(go.Scatter(x=distances, y=reference_amplitude, mode='lines', line=dict(dash='dash', color='black'),name=f'exp(−κd), κ={kappa:.3f}'))
    fig.update_layout(title=f'Exponential tail outside the {side} edge (E={E_n:.3f})', xaxis_title='Distance d from the well edge', yaxis_title='|φ|', yaxis_type='log', width=850, height=500)
    fig.show()

def plot_1d_time_evolution(phi_n, E_n, x_grid, N_t=200, t_max=10):
    #phi_n, E_n, x_grid as above, N_t: number of timesteps, t_max is max time

    X, T, Re_psi = real_time_evolved_1d(phi_n=phi_n, E_n=E_n, x_grid=x_grid, N_t=N_t, t_max=t_max)

    fig = go.Figure()
    fig.add_trace(go.Surface(x=X, y=T, z=Re_psi, colorscale='RdBu', cmid=0, colorbar=dict(title='Re(ψ)')))

    fig.update_layout(title=f'Real part of ψ(x,t) for state with E={E_n:.2f}', scene=dict(xaxis_title='x', yaxis_title='time t', zaxis_title='Re(ψ)'), width=900, height=600)

    fig.show()

def plot_superposition_prob(phi_1, E_1, phi_2, E_2, x_grid, t_max=10, N_t=300):
    #we take two eigenstates, and plot the superposition probability and plot as a function of x, 
    t_arr, prob = superposition_prob(phi_1=phi_1, phi_2=phi_2, E_1=E_1, E_2=E_2, x_grid=x_grid, t_max=t_max, N_t=N_t)

    fig = go.Figure(data=go.Surface(x=x_grid, y=t_arr, z=prob, colorscale='Viridis', colorbar=dict(title='|ψ|²')))
    fig.update_layout(title=f'Superposition of states E₁={E_1:.2f}, E₂={E_2:.2f}',scene=dict(xaxis_title='x',yaxis_title='time t',zaxis_title='|ψ|²'), width=900, height=600)

    fig.show()

def interactive_3d_well(a_init=4.0, b_init=3.0, c_init=5.0, V0_init=10.0, N_1d=75, N_3d=75):
    #solving to get energies, grids
    X3, Y3, Z3, energies, interps = solve_3d_well(a_init, b_init, c_init, V0_init, N_1d=N_1d, N_3d=N_3d)
    E_x, E_y, E_z = energies

    #defining sliders to be used on the plot
    a_s = widgets.FloatSlider(min=1.0, max=8.0, step=0.5, value=a_init, description='a')
    b_s = widgets.FloatSlider(min=1.0, max=8.0, step=0.5, value=b_init, description='b')
    c_s = widgets.FloatSlider(min=1.0, max=8.0, step=0.5, value=c_init, description='c')
    V0_s = widgets.FloatSlider(min=1.0, max=20.0, step=1.0, value=V0_init, description='V0')

    nx_s = widgets.IntSlider(min=0, max=len(E_x)-1, step=1, value=0, description='nx')
    ny_s = widgets.IntSlider(min=0, max=len(E_y)-1, step=1, value=0, description='ny')
    nz_s = widgets.IntSlider(min=0, max=len(E_z)-1, step=1, value=0, description='nz')
    core_iso_s = widgets.FloatSlider(min=0.1, max=0.9, step=0.1, value=0.5, description='core level') #high probability regions
    tail_iso_s = widgets.FloatLogSlider(value=0.03, base=10, min=-3, max=-0.7, step=0.1, description='tail level') #low probability regions

    out = widgets.Output()

    #keeping track of current geometry to avoid changing each time
    state = {'a': a_init, 'b': b_init, 'c': c_init, 'V0': V0_init, 'X3': X3, 'Y3': Y3, 'Z3': Z3, 'energies': energies, 'interps': interps}

    def update(a, b, c, V0, nx, ny, nz, core_iso_frac, tail_iso_frac):
        # we only recompute if things have changed
        if (a, b, c, V0) != (state['a'], state['b'], state['c'], state['V0']):
            X3_new, Y3_new, Z3_new, energies_new, interps_new = solve_3d_well(a, b, c, V0, N_1d=N_1d, N_3d=N_3d)
            state['a'], state['b'], state['c'], state['V0'] = a, b, c, V0
            state['X3'], state['Y3'], state['Z3'] = X3_new, Y3_new, Z3_new
            state['energies'] = energies_new
            state['interps'] = interps_new
            #updating max slider values
            E_x_new, E_y_new, E_z_new = energies_new
            nx_s.max = len(E_x_new)-1
            ny_s.max = len(E_y_new)-1
            nz_s.max = len(E_z_new)-1
        else:
            X3_new, Y3_new, Z3_new = state['X3'], state['Y3'], state['Z3']
            energies_new = state['energies']
            interps_new = state['interps']

        E_x_new, E_y_new, E_z_new = energies_new
        #restrict quantum numbers to new ranges
        nx = min(nx, len(E_x_new)-1)
        ny = min(ny, len(E_y_new)-1)
        nz = min(nz, len(E_z_new)-1)
        #re-solving equations using updated info
        phi_3d, prob_3d = compute_3d_state_seperable(X3_new, Y3_new, Z3_new, interps_new, (nx, ny, nz))
        E_total = E_x_new[nx] + E_y_new[ny] + E_z_new[nz]
        max_prob = prob_3d.max()
        core_iso_val = core_iso_frac * max_prob
        tail_iso_val = tail_iso_frac * max_prob

        #when integrating over the new grid, we may not have a total of 1, so we normalise the grid
        x_grid = X3_new[:, 0, 0]
        y_grid = Y3_new[0, :, 0]
        z_grid = Z3_new[0, 0, :]
        outside_box = ((np.abs(X3_new) > a / 2) | (np.abs(Y3_new) > b / 2) | (np.abs(Z3_new) > c / 2))
        total_probability = simpson(simpson(simpson(prob_3d, x=z_grid, axis=2), x=y_grid, axis=1),x=x_grid, axis=0)
        outside_probability = simpson(simpson(simpson(prob_3d * outside_box, x=z_grid, axis=2), x=y_grid, axis=1), x=x_grid, axis=0)
        p_outside = outside_probability / total_probability

        #we draw the twelve edges of the finite well so its obviously visually
        hx, hy, hz = a / 2, b / 2, c / 2
        vertices = [(-hx, -hy, -hz), (hx, -hy, -hz), (hx, hy, -hz), (-hx, hy, -hz), (-hx, -hy, hz),  (hx, -hy, hz),  (hx, hy, hz),  (-hx, hy, hz),]
        edges = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7)]
        edge_x, edge_y, edge_z = [], [], []
        for start, end in edges:
            for coordinate, values in enumerate((edge_x, edge_y, edge_z)):
                values.extend([vertices[start][coordinate], vertices[end][coordinate], None])

        with out:
            out.clear_output(wait=True) #light blue isosurface for low probs
            fig = go.Figure()
            fig.add_trace(go.Isosurface(
                x=X3_new.flatten(),
                y=Y3_new.flatten(),
                z=Z3_new.flatten(),
                value=prob_3d.flatten(),
                isomin=tail_iso_val,
                isomax=tail_iso_val,
                surface_count=1,
                colorscale='Blues',
                opacity=0.22,
                showscale=False,
                name='low-probability tail'
            ))
            fig.add_trace(go.Isosurface( #green isosurface for high probs
                x=X3_new.flatten(),
                y=Y3_new.flatten(),
                z=Z3_new.flatten(),
                value=prob_3d.flatten(),
                isomin=core_iso_val,
                isomax=core_iso_val,
                surface_count=1,
                colorscale='Viridis',
                opacity=0.82,
                showscale=False,
                name='high-probability core'
            ))
            fig.add_trace(go.Scatter3d( #dashed line to represent box
                x=edge_x, y=edge_y, z=edge_z, mode='lines',
                line=dict(color='black', width=5, dash='dash'),
                name='well boundary'
            ))
            fig.update_layout(#stats to show the quantum numbers, energy, prob outside of selected state
                title=(f'3D state ({nx},{ny},{nz})  E={E_total:.2f} 'f'P(outside)={p_outside:.2%}<br>' f'(a={a}, b={b}, c={c}, V₀={V0})'),
                scene=dict(xaxis_title='x', yaxis_title='y', zaxis_title='z'),
                width=700, height=600
            )
            fig.show()

    #connecting sliders to change when values change
    widgets.interactive_output(update,{'a': a_s, 'b': b_s, 'c': c_s, 'V0': V0_s, 'nx': nx_s, 'ny': ny_s, 'nz': nz_s, 'core_iso_frac': core_iso_s, 'tail_iso_frac': tail_iso_s})
    #4 sliders on top row, 5 on bottom
    ui = widgets.VBox([widgets.HBox([a_s, b_s, c_s, V0_s]), widgets.HBox([nx_s, ny_s, nz_s, core_iso_s, tail_iso_s]),out])
    display(ui)

def plot_3d_probability_outside_vs_energy(E_x, E_y, E_z, phi_x_int, phi_y_int, phi_z_int, X3, Y3, Z3, a, b, c, V0, max_states_per_axis=6):
    #plot 3d probability outside vs energy 
    nx_count = min(max_states_per_axis, len(E_x))
    ny_count = min(max_states_per_axis, len(E_y))
    nz_count = min(max_states_per_axis, len(E_z))

    #creating grids from the meshgrid input
    x_grid = X3[:, 0, 0]
    y_grid = Y3[0, :, 0]
    z_grid = Z3[0, 0, :]
    #creating phi's from interps input
    phi_x = np.column_stack([f(x_grid) for f in phi_x_int])
    phi_y = np.column_stack([f(y_grid) for f in phi_y_int])
    phi_z = np.column_stack([f(z_grid) for f in phi_z_int])

    nx, ny, nz = np.meshgrid(np.arange(nx_count), np.arange(ny_count), np.arange(nz_count), indexing='ij')
    #calculating energy levels and probabilities in each direction
    total_energy = (E_x[nx] + E_y[ny] + E_z[nz]).ravel()
    P_out_x = probability_outside_1d(x_grid, phi_x, E_x, a / 2)
    P_out_y = probability_outside_1d(y_grid, phi_y, E_y, b / 2)
    P_out_z = probability_outside_1d(z_grid, phi_z, E_z, c / 2)
    #total probabilty = 1, particle in box iff all components are in the box
    p_outside = (1 - ((1 - np.asarray(P_out_x)[nx]) * (1 - np.asarray(P_out_y)[ny]) * (1 - np.asarray(P_out_z)[nz]))).ravel()
    state_labels = np.array([f'({i},{j},{k})' for i, j, k in zip(nx.ravel(), ny.ravel(), nz.ravel())])
    bound = total_energy < V0
    #combining into 1 figure
    fig = go.Figure()
    for mask, name, colour, symbol in [(bound, 'Bound: E < V₀', '#168a4b', 'circle'), (~bound, 'Above V₀: E ≥ V₀', '#d1495b', 'x'),]: #green circle for bound states, red cross for unbound
        fig.add_trace(go.Scatter(
            x=total_energy[mask], y=p_outside[mask], mode='markers',
            customdata=state_labels[mask, np.newaxis],
            marker=dict(size=9, color=colour, symbol=symbol),
            hovertemplate=(
                'state = %{customdata[0]}<br>E = %{x:.3f}<br>'
                'P(outside) = %{y:.3%}<extra></extra>'
            ),
            name=name
        ))
    fig.add_vline(x=V0, line_dash='dash', line_color='black', annotation_text='V₀', annotation_position='top left') #dotted line to signify bound/unbound states
    fig.update_layout(
        title=f'3D exterior probability versus total energy (V₀ = {V0:g})',
        xaxis_title='Total energy Eₓ + Eᵧ + E_z',
        yaxis_title='P(outside 3D box)',
        yaxis_tickformat='.0%',
        width=850,
        height=550
    )
    fig.show()