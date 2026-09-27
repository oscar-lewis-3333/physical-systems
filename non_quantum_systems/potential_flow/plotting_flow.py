#static plots for various

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from model_flow import dipole_complex_potential, dipole_complex_velocity, uniform_complex_potential, vortex_complex_potential, vortex_complex_velocity


def plot_streamlines_and_equipotentials(x, y, phi, psi, *, u=None, v=None, stream_levels=20, potential_levels=20, stream_density=1.2, arrow_size=1.2,
    stagnation_points=None, show_equipotentials=True, title="Streamlines and equipotential lines", ax=None, show=True):
    #does what it says on the tin
    x = np.asarray(x)
    y = np.asarray(y)
    phi = np.ma.asarray(phi)
    psi = np.ma.asarray(psi)

    if phi.ndim != 2 or psi.ndim != 2:
        raise ValueError("phi and psi must be two-dimensional arrays.")
    if phi.shape != psi.shape:
        raise ValueError("phi and psi must have the same shape.")

    if (u is None) != (v is None):
        raise ValueError("Supply both velocity components u and v, or neither.")
    if u is not None:
        u = np.ma.asarray(u)
        v = np.ma.asarray(v)
        if u.shape != phi.shape or v.shape != phi.shape:
            raise ValueError("u and v must have the same shape as phi and psi.")
    if x.ndim == 1 and y.ndim == 1:
        if phi.shape != (y.size, x.size):
            raise ValueError(
                "For one-dimensional x and y, phi and psi must have shape (len(y), len(x)).")
    elif x.ndim == 2 and y.ndim == 2:
        if x.shape != phi.shape or y.shape != phi.shape:
            raise ValueError("Meshgrid arrays x and y must match phi and psi.")
    else:
        raise ValueError("x and y must both be one-dimensional or both be two-dimensional.")


    #apply one common mask so singular regions are omitted everywhere
    invalid = np.ma.getmaskarray(phi) | np.ma.getmaskarray(psi)
    invalid |= ~np.isfinite(np.ma.filled(phi, np.nan))
    invalid |= ~np.isfinite(np.ma.filled(psi, np.nan))

    if u is not None:
        invalid |= np.ma.getmaskarray(u) | np.ma.getmaskarray(v)
        invalid |= ~np.isfinite(np.ma.filled(u, np.nan))
        invalid |= ~np.isfinite(np.ma.filled(v, np.nan))
        u = np.ma.array(u, mask=invalid)
        v = np.ma.array(v, mask=invalid)

    phi = np.ma.array(phi, mask=invalid)
    psi = np.ma.array(psi, mask=invalid)

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))
    else:
        fig = ax.figure

    if u is None:
        ax.contour(x, y, psi, levels=stream_levels, colors="tab:blue", linestyles="-", linewidths=1.2)
    else:
        ax.streamplot(x, y, u, v, color="tab:blue", density=stream_density, linewidth=1.2, arrowsize=arrow_size)
    legend_handles = [Line2D([0], [0], color="tab:blue", linewidth=1.2, label=r"streamlines: $\psi=\mathrm{constant}$")]

    if show_equipotentials:
        ax.contour(x, y, phi, levels=potential_levels, colors="tab:orange", linestyles="--", linewidths=1.0)
        legend_handles.append(Line2D([0], [0], color="tab:orange", linestyle="--", linewidth=1.0, label=r"equipotential lines: $\phi=\mathrm{constant}$"))

    if stagnation_points is not None:
        points = np.asarray(stagnation_points)
        if np.iscomplexobj(points):
            points = np.atleast_1d(points)
            stagnation_x = points.real.ravel()
            stagnation_y = points.imag.ravel()
        else:
            points = np.asarray(stagnation_points, dtype=float)
            if points.ndim == 0:
                points = np.array([[points.item(), 0.0]])
            elif points.shape == (2,):
                points = points.reshape(1, 2)
            if points.ndim != 2 or points.shape[1] != 2:
                raise ValueError("stagnation_points must be a real or complex number, a sequence of complex numbers, or coordinate pairs with shape (n, 2).")
            stagnation_x = points[:, 0]
            stagnation_y = points[:, 1]

        point_label = "stagnation point" if stagnation_x.size == 1 else "stagnation points"
        point_handle = ax.scatter(stagnation_x, stagnation_y, s=60, color="black", edgecolors="white", linewidths=0.8, zorder=5, label=point_label)
        legend_handles.append(point_handle)

    ax.legend(handles=legend_handles)
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")
    ax.set_title(title)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(alpha=0.2)
    fig.tight_layout()

    if show:
        plt.show()

    return fig, ax


def plot_cylinder_flow(U=1.0, a=1.0, *, x_limits=(-4.0, 4.0), y_limits=(-4.0, 4.0), resolution=600, stream_density=1.2, arrow_size=1.2, ax=None, show=True):
    #plotting uniform potential flow around a circular cylinder
    if U <= 0:
        raise ValueError("U must be positive.")
    if a <= 0:
        raise ValueError("a must be positive.")
    if resolution < 2:
        raise ValueError("resolution must be at least 2.")

    D = 2 * np.pi * U * a**2
    x = np.linspace(*x_limits, resolution)
    y = np.linspace(*y_limits, resolution)
    X, Y = np.meshgrid(x, y)
    Z = X + 1j * Y

    #exclude interior from the fluid domain
    Z = np.ma.masked_where(np.abs(Z) <= a, Z)

    complex_potential = uniform_complex_potential(Z, U) + dipole_complex_potential(Z, D)
    complex_velocity = U + dipole_complex_velocity(Z, D)

    phi = complex_potential.real
    psi = complex_potential.imag
    u = complex_velocity.real
    v = -complex_velocity.imag

    fig, ax = plot_streamlines_and_equipotentials(X, Y, phi, psi, u=u, v=v, stream_density=stream_density, arrow_size=arrow_size, stagnation_points=[(-a, 0), (a, 0)],
        title=rf"Potential flow around a cylinder: $U={U}$, $a={a}$", ax=ax, show=False)

    cylinder = plt.Circle((0, 0), a, facecolor="0.92", edgecolor="black", linewidth=2.0, zorder=4)
    ax.add_patch(cylinder)

    legend = ax.get_legend()
    legend.set_loc("upper left")
    legend.set_bbox_to_anchor((1.02, 1.0))

    ax.set_xlim(*x_limits)
    ax.set_ylim(*y_limits)
    fig.tight_layout()

    if show:
        plt.show()

    return fig, ax


def plot_cylinder_pressure_distribution(U=1.0, a=1.0, Gamma=0.0, *, points=500, compare_zero_circulation=False, ax=None, show=True):
    #plot cylinders surface coefficient (with circulation option)
    if U <= 0:
        raise ValueError("U must be positive.")
    if a <= 0:
        raise ValueError("a must be positive.")
    if points < 2:
        raise ValueError("points must be at least 2.")

    theta = np.linspace(0, 2 * np.pi, points)
    circulation_term = Gamma / (2 * np.pi * a * U)
    pressure_coefficient = 1 - (-2 * np.sin(theta) + circulation_term) ** 2
    circulation_ratio = Gamma / (4 * np.pi * U * a)

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))
    else:
        fig = ax.figure

    ax.plot(theta, pressure_coefficient, color="tab:blue", linewidth=2.0, label=rf"$\Gamma/(4\pi Ua)={circulation_ratio:.2f}$")
    if compare_zero_circulation and not np.isclose(Gamma, 0.0):
        ax.plot(theta, 1 - 4 * np.sin(theta) ** 2, color="0.45", linewidth=1.5, linestyle="--", label=r"$\Gamma=0$")

    if np.abs(circulation_ratio) <= 1:
        first_angle = np.arcsin(circulation_ratio)
        stagnation_angles = np.mod([first_angle, np.pi - first_angle], 2 * np.pi,)
        if np.isclose(stagnation_angles[0], stagnation_angles[1]):
            stagnation_angles = stagnation_angles[:1]
        ax.scatter(stagnation_angles, np.ones_like(stagnation_angles), color="black", s=45, zorder=3, label="surface stagnation points")

    if np.isclose(Gamma, 0.0):
        ax.scatter([np.pi / 2, 3 * np.pi / 2], [-3, -3], color="tab:red", s=45, zorder=3, label="minimum pressure")
    ax.axhline(0, color="0.4", linewidth=1.0, linestyle="--")
    ax.set_xticks([0, np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi], [r"$0$", r"$\pi/2$", r"$\pi$", r"$3\pi/2$", r"$2\pi$"])
    ax.set_xlabel(r"Angular position, $\theta$")
    ax.set_ylabel(r"Pressure coefficient, $C_p$")
    ax.set_title("Pressure distribution around the cylinder surface")
    ax.set_xlim(0, 2 * np.pi)
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()

    if show:
        plt.show()

    return fig, ax


def plot_circulating_cylinder_flow(U=1.0, a=1.0, Gamma=0.0, *, x_limits=(-4.0, 4.0), y_limits=(-4.0, 4.0), resolution=600, stream_density=1.4, arrow_size=1.2, ax=None, show=True):
    #plotting cylinder flow with a given circulation
    if U <= 0:
        raise ValueError("U must be positive.")
    if a <= 0:
        raise ValueError("a must be positive.")
    if resolution < 2:
        raise ValueError("resolution must be at least 2.")

    D = 2 * np.pi * U * a**2
    x = np.linspace(*x_limits, resolution)
    y = np.linspace(*y_limits, resolution)
    X, Y = np.meshgrid(x, y)
    Z = X + 1j * Y
    Z = np.ma.masked_where(np.abs(Z) <= a, Z)

    complex_potential = uniform_complex_potential(Z, U) + dipole_complex_potential(Z, D) + vortex_complex_potential(Z, Gamma)
    complex_velocity = U + dipole_complex_velocity(Z, D) + vortex_complex_velocity(Z, Gamma)

    phi = complex_potential.real
    psi = complex_potential.imag
    u = complex_velocity.real
    v = -complex_velocity.imag

    #roots of U z^2 - i Gamma z/(2 pi) - U a^2 = 0. 
    roots = np.roots([U, -1j * Gamma / (2 * np.pi), -U * a**2])
    stagnation_points = []
    for root in roots:
        if np.abs(root) >= a * (1 - 1e-8):
            if not any(np.isclose(root, point) for point in stagnation_points):
                stagnation_points.append(root)

    fig, ax = plot_streamlines_and_equipotentials(X, Y, phi, psi, u=u, v=v, stream_density=stream_density, arrow_size=arrow_size,
        stagnation_points=np.asarray(stagnation_points, dtype=complex), show_equipotentials=False,
        title=rf"Cylinder flow with circulation: $U={U}$, $a={a}$, $\Gamma={Gamma:.3g}$",
        ax=ax, show=False)

    cylinder = plt.Circle((0, 0), a, facecolor="0.92", edgecolor="black", linewidth=2.0, zorder=4)
    ax.add_patch(cylinder)

    legend = ax.get_legend()
    legend.set_loc("upper left")
    legend.set_bbox_to_anchor((1.02, 1.0))

    ax.set_xlim(*x_limits)
    ax.set_ylim(*y_limits)
    fig.tight_layout()

    if show:
        plt.show()

    return fig, ax
