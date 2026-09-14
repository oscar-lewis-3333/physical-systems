# Bead on a Rotating Hoop

![Animation of a bead on a rotating hoop](images/bead_on_rotating_hoop.gif)

This project investigates the mechanics behind the above GIF from a numerical and analytic perspective. Start by deriving the equations of motion using Euler-Lagrange, solving numerically, before discussing equilibrium points and their behaviour, and the effective energy and potential of the system. Ends with the above render. 

## Assumptions

The model assumes:

- bead is a point-mass
- no forces acting between bead and hoop (friction etc.)
- hoop is rigid and spherical
- constant angular velocity
- no damping

For the angular velocity to stay constant, must be some form of external motor to maintain this, hence the energy may not be conserved. The effective energy is the quantity which is conserved in time in this scenario.

## Model

The hoop has radius \(R\) and rotates about its vertical axis with constant angular velocity \(\Omega\). The bead’s position is described by the angle \(\theta\), measured from the downward vertical (opposite to convention). The coordinates describing the bead are:

\[
x=R\sin\theta\cos(\Omega t), \qquad
y=R\sin\theta\sin(\Omega t), \qquad
z=-R\cos\theta.
\]

The resulting equation of motion (using Euler-Lagrange) is:

\[
\ddot{\theta}
=
\sin\theta
\left(
\Omega^2\cos\theta-\frac{g}{R}
\right).
\]

Because this equation is a nonlinear ODE, its trajectories are calculated numerically using SciPy's `solve_ivp`.

## Critical rotation speed

The system has the critical angular velocity:

\[
\Omega_c=\sqrt{\frac{g}{R}}.
\]

The stability of equilibrium points are dependent on the ratio:

\[
\lambda=\frac{\Omega}{\Omega_c}.
\]

- For \(\lambda<1\), the bottom equilibrium \(\theta=0\) is stable.
- At \(\lambda=1\), the bottom equilibrium becomes degenerate.
- For \(\lambda>1\), the bottom equilibrium becomes unstable and two stable equilibria appear:

\[
\theta_e
=
\pm\cos^{-1}\left(\frac{1}{\lambda^2}\right).
\]

This bifurcation is plotted in the notebook. The equilibrium at \(\theta=\pi\) is not plotted, but remains unstable for every value of \(\Omega\).

## Effective energy

The conserved reduced energy is:

\[
E_{\mathrm{eff}}
=
\frac12mR^2\dot{\theta}^2
-\frac12mR^2\Omega^2\sin^2\theta
-mgR\cos\theta.
\]

The effective potential is therefore:

\[
U_{\mathrm{eff}}(\theta)
=
-\frac12mR^2\Omega^2\sin^2\theta
-mgR\cos\theta+C.
\]

The project compares the numerical energy against its initial value to check the accuracy of the integration. It also uses the effective potential to identify allowed regions, equilibrium points and turning points.

## Small oscillations

Linearising around the stable equilibrium branches gives the dimensionless frequencies

\[
\frac{\omega_-}{\Omega_c}
=
\sqrt{1-\lambda^2},
\qquad \lambda<1,
\]

and

\[
\frac{\omega_+}{\Omega_c}
=
\sqrt{\lambda^2-\lambda^{-2}},
\qquad \lambda>1.
\]

Both frequencies approach zero at the bifurcation, producing a critical slowing down effect which is visible on the trajectory plot at(\lambda=1\)

## Project structure

```text
bead_on_rotating_hoop/
├── derivation.ipynb
├── solution_bead.py
├── plotting_bead.py
├── model.py
├── render.py
├── media/
│   └── bead_on_rotating_hoop.gif
│   └── bead_on_rotating_hoop.mp4
└── README.md
```

- `derivation.ipynb` — deriving equations of motion, interpreting plots/solution
- `solution_bead.py` — solving equation of motion, effective potential/energy functions
- `plotting_bead.py` — all static plots used in the notebook
- `model.py` — geometry of hoop/bead in cartesian co-ordinates
- `render.py` — 3D render of the system
