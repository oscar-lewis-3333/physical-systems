from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FFMpegWriter, FuncAnimation, PillowWriter
from matplotlib.patches import Circle, FancyArrowPatch


def render_lifting_cylinder(U, a, Gamma, rho, cylinder_mass_per_length, *, duration=3.0, fps=30, initial_height=0.0, initial_vertical_speed=0.0, angular_speed=None, initial_rotation_angle=0.0, include_added_mass=True, output_path=None, show=True):
    #animate a cylinder moving under the lift effect (Kutta-Joukowski)
    if U <= 0:
        raise ValueError("U must be positive.")
    if a <= 0:
        raise ValueError("a must be positive.")
    if rho <= 0:
        raise ValueError("rho must be positive.")
    if cylinder_mass_per_length <= 0:
        raise ValueError("cylinder_mass_per_length must be positive.")
    if duration <= 0:
        raise ValueError("duration must be positive.")
    if fps <= 0:
        raise ValueError("fps must be positive.")

    lift_per_length = -rho * U * Gamma
    added_mass_per_length = rho * np.pi * a**2 if include_added_mass else 0.0
    effective_mass_per_length = cylinder_mass_per_length + added_mass_per_length
    vertical_acceleration = lift_per_length / effective_mass_per_length
    if angular_speed is None:
        angular_speed = Gamma / (2 * np.pi * a**2)

    frame_count = max(2, int(np.ceil(duration * fps)) + 1)
    times = np.linspace(0.0, duration, frame_count)
    heights = initial_height + initial_vertical_speed * times + 0.5 * vertical_acceleration * times**2

    vertical_margin = 1.8 * a
    y_limits = (min(initial_height, heights.min()) - vertical_margin, max(initial_height, heights.max()) + vertical_margin,)
    x_limits = (-3.5 * a, 3.5 * a)

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.set_xlim(*x_limits)
    ax.set_ylim(*y_limits)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")
    ax.set_title("Cylinder motion under circulation-induced lift")
    ax.grid(alpha=0.2)

    arrow_x = np.linspace(x_limits[0] + 0.4 * a, x_limits[1] - 0.4 * a, 8)
    arrow_y = np.linspace(y_limits[0] + 0.4 * a, y_limits[1] - 0.4 * a, 8)
    X_arrow, Y_arrow = np.meshgrid(arrow_x, arrow_y)
    ax.quiver(X_arrow, Y_arrow, np.ones_like(X_arrow), np.zeros_like(Y_arrow), color="0.75", alpha=0.65, angles="xy", scale_units="xy", scale=1.8 / a, width=0.003, zorder=0)

    trail, = ax.plot([], [], color="tab:purple", linewidth=1.5, alpha=0.75)
    cylinder = Circle((0.0, initial_height), a, facecolor="0.90", edgecolor="black", linewidth=2.0, zorder=3)
    ax.add_patch(cylinder)

    lift_direction = np.sign(lift_per_length)
    lift_arrow_length = 1.2 * a
    lift_arrow = FancyArrowPatch((0.0, initial_height), (0.0, initial_height + lift_direction * lift_arrow_length), arrowstyle="-|>", mutation_scale=18, color="tab:red",
        linewidth=2.0, zorder=4)
    ax.add_patch(lift_arrow)

    #tangent arrow to travel around boundary
    surface_arrow_length = 0.65 * a

    def rotation_geometry(time, height):
        angle = initial_rotation_angle + angular_speed * time
        surface_point = np.array([a * np.cos(angle), height + a * np.sin(angle)])
        direction = np.sign(angular_speed) * np.array([-np.sin(angle), np.cos(angle)])
        start = surface_point - 0.45 * surface_arrow_length * direction
        end = surface_point + 0.55 * surface_arrow_length * direction
        return start, end, surface_point

    rotation_start, rotation_end, rotation_point = rotation_geometry(0.0, initial_height)
    rotation_arrow = FancyArrowPatch(rotation_start, rotation_end, arrowstyle="-|>", mutation_scale=16, color="tab:green", linewidth=2.0, zorder=5)
    ax.add_patch(rotation_arrow)
    rotation_marker, = ax.plot([rotation_point[0]], [rotation_point[1]], marker="o", markersize=5, color="tab:green", zorder=6)
    time_text = ax.text(0.02, 0.96, "", transform=ax.transAxes, va="top")
    parameter_text = ax.text(0.02, 0.88,"\n".join([rf"$\Gamma={Gamma:.3g}$", rf"$L={lift_per_length:.3g}$", rf"$a_y={vertical_acceleration:.3g}$",
                    rf"$\Omega={angular_speed:.3g}$"]), transform=ax.transAxes, va="top")

    if np.isclose(lift_per_length, 0.0):
        lift_arrow.set_visible(False)
    if np.isclose(angular_speed, 0.0):
        rotation_arrow.set_visible(False)
        rotation_marker.set_visible(False)

    def update(frame):
        height = heights[frame]
        cylinder.center = (0.0, height)
        trail.set_data(np.zeros(frame + 1), heights[: frame + 1])
        lift_arrow.set_positions((0.0, height), (0.0, height + lift_direction * lift_arrow_length))
        rotation_start, rotation_end, rotation_point = rotation_geometry(times[frame], height)
        rotation_arrow.set_positions(rotation_start, rotation_end)
        rotation_marker.set_data([rotation_point[0]], [rotation_point[1]])
        time_text.set_text(rf"$t={times[frame]:.2f}$")
        return (cylinder, trail, lift_arrow, rotation_arrow, rotation_marker, time_text, parameter_text)

    animation = FuncAnimation(fig, update, frames=frame_count, interval=1000 / fps, blit=True)
    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        suffix = output_path.suffix.lower()
        if suffix == ".gif":
            writer = PillowWriter(fps=fps)
        elif suffix == ".mp4":
            writer = FFMpegWriter(fps=fps, bitrate=2400)
        else:
            raise ValueError("output_path must end in .gif or .mp4.")
        animation.save(output_path, writer=writer)

    if show:
        plt.show()

    return animation, fig, ax
