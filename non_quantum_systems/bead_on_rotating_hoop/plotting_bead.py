import matplotlib.pyplot as plt
import numpy as np
from solution_bead import dimensionless_potential, dimensionless_inital_energy

#file only for making nice looking plots. there is nothing special happening here mathematically
def plot_sol(sol):
    theta = sol.y[0]
    theta_dot = sol.y[1]

    fig, axes = plt.subplots(2, 1, figsize=(12, 8))
    #plotting of theta over time
    axes[0].plot(sol.t, theta)
    axes[0].set_ylabel(r"$\theta(t)$")
    axes[0].grid()
    #plotting of theta dot over time
    axes[1].plot(sol.t, theta_dot)
    axes[1].set_xlabel('Time')
    axes[1].set_ylabel(r"$ \dot {\theta}(t)$")
    axes[1].grid()

    plt.tight_layout()
    plt.show()

    return fig

def phase_plane_comparison(sol_1, sol_2, sol_3, Omega_values, R, g):
    solutions = [sol_1, sol_2, sol_3]
    Omega_c = np.sqrt(g / R)
    if len(Omega_values) != len(solutions):
        raise ValueError("Supply one angular velocity for each solution.")

    theta_limit = 1.5 * max(np.max(np.abs(sol.y[0])) for sol in solutions)
    theta_dot_limit = 1.5 * max(np.max(np.abs(sol.y[1])) for sol in solutions)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), sharex=True, sharey=True)

    for ax, sol, Omega in zip(axes, solutions, Omega_values):
        theta = sol.y[0]
        theta_dot = sol.y[1]
        ratio = abs(Omega) / Omega_c

        ax.plot(theta, theta_dot, linewidth=1.5, label="trajectory")
        ax.scatter(theta[0], theta_dot[0], color="black", marker="*", s=70, zorder=4, label="initial state")

        if np.isclose(ratio, 1.0):
            regime = "critical"
            ax.scatter(0, 0, color="tab:orange", marker="D", s=45, zorder=5, label="critical equilibrium")
        elif ratio < 1.0:
            regime = "subcritical"
            ax.scatter(0, 0, color="tab:green", marker="o", s=45, zorder=5, label="stable equilibrium")
        else:
            regime = "supercritical"
            theta_e = np.arccos(g / (R * Omega**2))
            ax.scatter(0, 0, color="tab:red", marker="x", s=55, zorder=5, label="unstable equilibrium")
            ax.scatter([-theta_e, theta_e], [0, 0], color="tab:green", marker="o", s=45, zorder=5, label=r"stable equilibria $\pm\theta_e$")

        ax.axhline(0, color="0.75", linewidth=0.8)
        ax.axvline(0, color="0.75", linewidth=0.8)
        ax.set_title(rf"$\Omega/\Omega_c={ratio:.2f}$" + f" ({regime})")
        ax.set_xlabel(r"$\theta$ (rad)")
        ax.grid(alpha=0.3)

    axes[0].set_ylabel(r"$\dot{\theta}$ (rad s$^{-1}$)")
    axes[0].set_xlim(-theta_limit, theta_limit)
    axes[0].set_ylim(-theta_dot_limit, theta_dot_limit)
    fig.suptitle("Phase-plane comparison across the critical angular velocity")

    #trajectories often overlap most of the 'screen', so legend always looked awkwardly placed. adjust to just place one legend at the bottom
    legend_entries = {}
    for ax in axes:
        handles, labels = ax.get_legend_handles_labels()
        for handle, label in zip(handles, labels):
            legend_entries.setdefault(label, handle)

    fig.legend(legend_entries.values(), legend_entries.keys(), loc="lower center", bbox_to_anchor=(0.5, 0.01), ncol=3, fontsize=8)

    fig.tight_layout(rect=(0, 0.17, 1, 0.93))
    plt.show()

    return fig

def plot_energy(sol, m, R, g, Omega):

    from solution_bead import energy

    theta = sol.y[0]
    theta_dot = sol.y[1]
    energy_values = energy(y=[theta, theta_dot], m=m, R=R, g=g, Omega=Omega)
    E0 = energy_values[0]

    energy_error = (energy_values - E0) / (m * g * R)

    margin = 0.05 * max(abs(E0), m * g * R)

    fig, axes = plt.subplots(2, 1, figsize=(12, 8))
    #energy over time
    axes[0].plot(sol.t, energy_values)
    axes[0].set_ylabel('Energy')
    axes[0].set_ylim(E0 - margin, E0 + margin)
    axes[0].set_title('Effective energy in time')
    axes[0].grid()

    #difference in energy from its original value over time

    axes[1].plot(sol.t, energy_error)
    axes[1].set_ylabel(r'$(E-E(0))/(mgR)$')
    axes[1].set_xlabel('Time')
    axes[1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    axes[1].text(0.02, 0.95, rf"max error = ${np.max(np.abs(energy_error)):.2e}$", transform=axes[1].transAxes, va='top')
    axes[1].grid()

    plt.tight_layout()
    plt.show()

    return fig

def plot_dimensionless_potential(ratios, Omega_c):
    ratios = np.asarray(ratios, dtype=float)
    if ratios.size == 0:
        raise ValueError("Supply at least one value of Omega/Omega_c.")

    theta = np.linspace(-np.pi, np.pi, 500)
    fig, axes = plt.subplots(1, len(ratios), figsize=(4.7 * len(ratios), 4.8), sharex=True, sharey=True, squeeze=False)
    axes = axes.ravel()

    for ax, ratio in zip(axes, ratios):
        Omega = ratio * Omega_c
        U = dimensionless_potential(theta, Omega, Omega_c)
        speed_ratio = abs(ratio)

        ax.plot(theta, U, linewidth=2, label=r"$u(\theta)$")
        ax.axhline(0, color="0.65", linewidth=0.8)
        ax.axvline(0, color="0.65", linewidth=0.8)

        #theta = /pm pi represents the same unstable equilibrium at the top.
        U_top = dimensionless_potential(np.pi, Omega, Omega_c)
        ax.scatter([-np.pi, np.pi], [U_top, U_top], color="tab:red", marker="x", s=55, zorder=5, label=r"unstable top, $\theta=\pm\pi$",)

        if np.isclose(speed_ratio, 1.0):
            regime = "critical"
            ax.scatter(0, 0, color="tab:orange", marker="D", s=45, zorder=5, label=r"critical point, $\theta=0$")
        elif speed_ratio < 1.0:
            regime = "subcritical"
            ax.scatter(0, 0, color="tab:green", marker="o", s=45, zorder=5, label=r"stable bottom, $\theta=0$")
        else:
            regime = "supercritical"
            theta_e = np.arccos(1 / speed_ratio**2)
            U_equilibrium = dimensionless_potential(theta_e, Omega, Omega_c)
            ax.scatter(0, 0, color="tab:red", marker="x", s=55, zorder=5, label=r"unstable bottom, $\theta=0$")
            ax.scatter([-theta_e, theta_e], [U_equilibrium, U_equilibrium], color="tab:green", marker="o", s=45, zorder=5, label=r"stable points, $\theta=\pm\theta_e$")

        ax.set_title(rf"$|\Omega|/\Omega_c={speed_ratio:.2f}$" + f" ({regime})")
        ax.set_xlabel(r"$\theta$ (rad)")
        ax.set_xticks([-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi], [r"$-\pi$", r"$-\pi/2$", r"$0$", r"$\pi/2$", r"$\pi$"])
        ax.set_xlim(-np.pi - 0.08, np.pi + 0.08)
        ax.grid(alpha=0.3)

    axes[0].set_ylabel(r"$u(\theta)=[U_{\mathrm{eff}}(\theta)-U_{\mathrm{eff}}(0)]/(mgR)$")
    fig.suptitle("Dimensionless effective-potential landscape")

    #'same' legend at the bottom as for phase planes
    legend_entries = {}
    for ax in axes:
        handles, labels = ax.get_legend_handles_labels()
        for handle, label in zip(handles, labels):
            legend_entries.setdefault(label, handle)

    fig.legend(legend_entries.values(), legend_entries.keys(), loc="lower center", bbox_to_anchor=(0.5, 0.01), ncol=3, fontsize=8)

    fig.tight_layout(rect=(0, 0.16, 1, 0.92))
    plt.show()
    return fig

def plot_dimensionless_energy(y0, Omega, Omega_c):
    energy_0 = dimensionless_inital_energy(y0, Omega, Omega_c)
    theta_0 = y0[0]
    theta_dot_0 = y0[1]
    ratio = abs(Omega / Omega_c)

    theta = np.linspace(-np.pi, np.pi, 2000)
    U = dimensionless_potential(theta, Omega, Omega_c)
    U_initial = dimensionless_potential(theta_0, Omega, Omega_c)

    #locating intersections u(\theta)=e_0 by linear interpolation.
    difference = U - energy_0
    crossing_indices = np.where(np.signbit(difference[:-1]) != np.signbit(difference[1:]))[0]
    turning_points = []
    for index in crossing_indices:
        theta_left, theta_right = theta[index], theta[index + 1]
        difference_left, difference_right = difference[index], difference[index + 1]
        root = theta_left - difference_left * (theta_right - theta_left) / (difference_right - difference_left)
        turning_points.append(root)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(theta, U, linewidth=2, label=r"effective potential $u(\theta)$")
    ax.axhline(energy_0, color="tab:red", linestyle="--", linewidth=1.5, label=rf"initial energy $e_0={energy_0:.3f}$")
    #motion can only classically exist when e - u(\theta) >=0
    ax.fill_between(theta, U, energy_0, where=energy_0 >= U, color="tab:green", alpha=0.12, interpolate=True, label=r"classically allowed: $e_0\geq u(\theta)$")

    ax.scatter(theta_0, U_initial, color="black", marker="*", s=85, zorder=6, label=rf"initial position $\theta_0={theta_0:.2f}$")

    if not np.isclose(theta_dot_0, 0.0):
        ax.vlines(theta_0, U_initial, energy_0, color="black", linestyle=":", linewidth=1.2, label="initial kinetic contribution")

    if turning_points:
        ax.scatter(turning_points, np.full(len(turning_points), energy_0), color="tab:purple", marker="o", s=38, zorder=5, label=r"turning points: $e_0=u(\theta)$")

    U_top = dimensionless_potential(np.pi, Omega, Omega_c)
    ax.scatter([-np.pi, np.pi], [U_top, U_top], color="tab:red", marker="x", s=55, zorder=5, label=r"unstable top, $\theta=\pm\pi$")

    if np.isclose(ratio, 1.0):
        regime = "critical"
        ax.scatter(0, 0, color="tab:orange", marker="D", s=45, zorder=5, label=r"critical point, $\theta=0$")
    elif ratio < 1.0:
        regime = "subcritical"
        ax.scatter(0, 0, color="tab:green", marker="o", s=45, zorder=5, label=r"stable bottom, $\theta=0$")
    else:
        regime = "supercritical"
        theta_e = np.arccos(1 / ratio**2)
        U_equilibrium = dimensionless_potential(theta_e, Omega, Omega_c)
        ax.scatter(0, 0, color="tab:red", marker="x", s=55, zorder=5, label=r"unstable bottom, $\theta=0$")
        ax.scatter([-theta_e, theta_e], [U_equilibrium, U_equilibrium], color="tab:green", marker="o", s=45, zorder=5, label=r"stable points, $\theta=\pm\theta_e$")

    ax.axhline(0, color="0.65", linewidth=0.8)
    ax.axvline(0, color="0.65", linewidth=0.8)
    ax.set_xlabel(r"$\theta$ (rad)")
    ax.set_ylabel("Dimensionless effective energy")
    ax.set_title("Allowed motion in the effective potential\n" + rf"$|\Omega|/\Omega_c={ratio:.2f}$ ({regime})")
    ax.set_xticks([-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi], [r"$-\pi$", r"$-\pi/2$", r"$0$", r"$\pi/2$", r"$\pi$"])
    ax.set_xlim(-np.pi - 0.08, np.pi + 0.08)
    ax.grid(alpha=0.3)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=3, fontsize=8)
    fig.tight_layout()
    plt.show()

    return fig

    
def plot_bifurcation(ratio_max=2.5, number_of_points=1000):
    if ratio_max <= 1:
        raise ValueError("ratio_max must be greater than 1.")

    #critcal point at \lambda=1, so ratio above and below this point. both of these are \lambda just defined on different domains
    ratio_below = np.linspace(0, 1, number_of_points)
    ratio_above = np.linspace(1, ratio_max, number_of_points)

    #non-zero equilibria which only exist for ratio >= 1
    theta_e = np.arccos(1 / ratio_above**2)

    fig, ax = plt.subplots(figsize=(8, 6))

    #theta=0 is stable below critical point
    ax.plot(ratio_below, np.zeros_like(ratio_below), color="tab:green", linestyle="-", linewidth=2, label=r"stable $\theta=0$")

    #theta=0 is unstable past the critical point
    ax.plot(ratio_above, np.zeros_like(ratio_above), color="tab:red", linestyle="--", linewidth=2, label=r"unstable $\theta=0$")

    #plotting \pm \theta_{e}
    ax.plot(ratio_above, theta_e, color="tab:green", linestyle="-", linewidth=2, label=r"stable $\theta=+\theta_e$")
    ax.plot(ratio_above, -theta_e, color="tab:green", linestyle="-", linewidth=2, label=r"stable $\theta=-\theta_e$")
    #critical point at \lambda=1
    ax.scatter(1, 0, color="tab:orange", marker="D", s=65, zorder=5, label=r"critical point $\lambda=1$")
    ax.axvline(1, color="0.5", linestyle=":", linewidth=1)
    ax.set_xlabel(r"$\lambda=|\Omega|/\Omega_c$")
    ax.set_ylabel(r"Equilibrium angle $\theta_e$ (rad)")
    ax.set_title("Pitchfork bifurcation of the rotating-hoop equilibria")
    ax.grid(alpha=0.3)
    ax.legend()

    fig.tight_layout()
    plt.show()

    return fig


    
