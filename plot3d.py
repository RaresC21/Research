import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from Polytope import Polytope
from symmetric import *

def plot_plane(plane, plt):
    point = plane.pt_on_plane
    normal = plane.perp_vec

    # a plane is a*x+b*y+c*z+d=0
    # [a,b,c] is the normal. Thus, we have to calculate
    # d and we're set
    point = np.array(point)
    d = np.dot(point * -1, normal)

    # create x,y
    xx, yy = np.meshgrid(range(-2,3), range(-2,3))

    # calculate corresponding z
    z = (-normal[0] * xx - normal[1] * yy - d) * 1. /normal[2]

    # plot the surface
    plt.plot_surface(xx, yy, z)

def plot_polytope(polytope, plt):
    for p in polytope.planes:
        plot_plane(p, plt)

plt3d = plt.figure().gca(projection='3d')
poly = symmetric(3)
plot_polytope(poly, plt3d)
plt.show()
