import numpy as np

#nothing special is happening here, just setting up our model for the 3d render
def hoop_coordinates(alpha, t, R, Omega):
    #all of this is spherical coordinates
    phi = Omega * t 

    x = R * np.sin(alpha) * np.cos(phi)
    y = R * np.sin(alpha) * np.sin(phi)
    z = -R * np.cos(alpha)

    return x, y, z


def bead_coordinates(theta, t, R, Omega):
    return hoop_coordinates(theta, t, R, Omega)