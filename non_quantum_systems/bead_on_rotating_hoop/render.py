import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from model import hoop_coordinates, bead_coordinates

#there is nothing cool mathematically going on here, but just making a 3d render look better

def animate_solution(sol, R, Omega, fps=30, frame_step=2, trail_length=40):
    times = sol.t[::frame_step]
    theta = sol.y[0, ::frame_step]

    alpha = np.linspace(-np.pi, np.pi, 300)
    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(111, projection="3d")

    limit = 1.15 * R

    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_zlim(-limit, limit)
    ax.set_box_aspect((1, 1, 1))

    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    ax.set_zlabel("$z$")
    ax.view_init(elev=20, azim=35)

    #fix the vertical rotation axis
    ax.plot([0, 0], [0, 0], [-limit, limit], linestyle="--", color="grey", label="Rotation axis")
    hoop_line, = ax.plot([], [], [], color="royalblue", linewidth=2, label="Hoop")
    trail_line, = ax.plot([], [], [], color="crimson", alpha=0.4)
    bead_point, = ax.plot([], [], [], "o", color="crimson", markersize=9, label="Bead")

    title = ax.set_title("")
    ax.legend(loc="upper left")

    def update(frame):
        time = times[frame]

        #update hoop
        x_hoop, y_hoop, z_hoop = hoop_coordinates(alpha, time, R, Omega)
        hoop_line.set_data_3d(x_hoop, y_hoop, z_hoop)

        #update bead
        x_bead, y_bead, z_bead = bead_coordinates(theta[frame], time, R, Omega)
        bead_point.set_data_3d([x_bead], [y_bead], [z_bead])

        #show most recent section of path
        start = max(0, frame - trail_length)

        x_trail, y_trail, z_trail = bead_coordinates(theta[start:frame + 1], times[start:frame + 1], R, Omega)
        trail_line.set_data_3d(x_trail, y_trail, z_trail)

        title.set_text(f"t = {time:.2f} s, theta = {theta[frame]:.2f} rad")

        return hoop_line, bead_point, trail_line, title

    animation = FuncAnimation(fig, update, frames=len(times), interval=1000 / fps, blit=False)

    return fig, animation