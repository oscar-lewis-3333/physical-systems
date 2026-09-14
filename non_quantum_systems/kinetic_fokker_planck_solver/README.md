# Kinetic Fokker–Planck Solver

A finite-difference solver for the following 1D kinetic Fokker–Planck equation in phase space:

\[
\frac{\partial f}{\partial t}
=-v\frac{\partial f}{\partial x}
-F(x)\frac{\partial f}{\partial v}
+C\frac{\partial^2 f}{\partial v^2}.
\]

The code advances the distribution `f(x, v, t)` using operator splitting at each time step, to then solve implicitly in both x, v directions  Both stages use upwind finite differences and LAPACKE for solving resulting linear systems.

## Files

- `specification.pdf` - the question
- `kinetic_fokker_planck_solver.c` — main solver.
- `input.txt` — grid, time and diffusion parameters for the example run.
- `coefficients.txt` — given values of the force function `F(x)`.
- `output.txt` — solution based on input, coefficients
- `Makefile` — build, run and clean commands
- `.vscode/tasks.json` — VSCode build and run tasks

## Input format

`input.txt` contains seven values, one per line:

1. `N_x` — number of spatial grid points.
2. `N_v` — number of velocity grid points.
3. `t_f` — final simulation time.
4. `L` — length of the spatial domain `[0, L]`.
5. `v_m` — velocity bound defining `[-v_m, v_m]`.
6. `C` — velocity-diffusion coefficient.
7. `I_min` — minimum number of time steps.

`coefficients.txt` must then contain `N_x` sampled values of `F(x)`, else a safe exit error. The solver writes `N_x × N_v` values to `output.txt`, ordered by the spatial variable then the velocity variable

## Numerical notes
Time-step chosen to respect the CFL condition before each direction is solved implicitly. Boundary conditions enforced directly inside the banded systems, due to applying at any other time causing issues closer to boundary (for example, the boundary conditions may be enforced but may be a jump near the boundary). It is not a general use PDE package, moreso a showcase of operator splitting and a solution using finite-difference (linear algebra) methods.
