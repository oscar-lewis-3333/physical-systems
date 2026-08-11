# Quantum Tunnelling — Finite Potential Well
A numerical and visual exploration of quantum tunnelling using the problem of the finite potential well alongside time-independent and time-dependent Schrödinger equations.

The project starts with the one-dimensional time-independent case, before discussing the 1D time-dependent case, then the 3D case. We use eigenvalue methods to calculate bound states, eigenstates, and then the energy spectra and probability density across the problem.

## Features

- Numerical solution of the 1D time-independent Schrödinger equation
- Visualisation of the quantised bound-states and corresponding eigenstates
- Calculation of probability densities and the probability of an eigenstate existing outside the well
- Comparison of numerical wavefunction tails with the theoretical exponential decay
  $$
  |\phi(x)| \propto e^{-\kappa d},
  \qquad
  \kappa = \frac{\sqrt{2m(V_0-E)}}{\hbar}
  $$
- Time evolution of stationary states and two-state superpositions, alongside the probability distrubution as a function of position and time
- Interactive 3D probability-density isosurfaces with adjustable quantum numbers
- 3D energy versus exterior-probability plot, distinguishing states below and above the chosen $V_0$ threshold

## Project structure

```text
physical_systems/
└── Non-Quantum
└── Quantum
    └── quantum_tunnelling/
        ├── finite_potential_well.ipynb  # Main analysis and visualisations
        ├── solving.py                   # Numerical solvers and wavefunction utilities
        └── plotting.py                  # Plotting and interactive visualisation functions
```

## Model

For the 1D finite well,

$$
V(x)=
\begin{cases}
0, & |x|<a, \\
V_0, & |x|\geq a.
\end{cases}
$$

Bound states satisfy $E<V_0$. In classical mechanics, the particle cannot exist outside the well, since it would have Kinetic Energy $T = E - V_0 < 0$, but we see that in quantum mechanics we have an extension into the forbidden region with exponential decay.

For the 3D numerical extension, the project uses a separable potential:

$$
V(x,y,z)=V_x(x)+V_y(y)+V_z(z).
$$

This gives product-state solutions:

$$
\phi_{n_x,n_y,n_z}(x,y,z)=\phi_{n_x}(x)\phi_{n_y}(y)\phi_{n_z}(z),
$$

with total energy

$$
E_{n_x,n_y,n_z}=E_{n_x}^{(x)}+E_{n_y}^{(y)}+E_{n_z}^{(z)}.
$$

## Key observations

- Quantised energy levels observed
- Higher-energy bound states generally have greater probability of existing outside the well
- In the classically forbidden region, the wavefunction has exponential decay with respect to spatial variable
- A stationary energy eigenstate has a time-independent probability density.
- A superposition of energy eigenstates has an sinusiodal probability distrubution when restricted to time
- In 3D, nodal planes divide a single continuous wavefunction into various regions of high-probability of existence

## Limitations

- The calculations use dimensionless units, with $\hbar=m=1$.
- The numerical domain is finite, so the exponentially decaying tails are approximated on a finite grid.
- The 3D solver uses a separable approximation rather than solving the full non-separable rectangular finite-well problem directly.
