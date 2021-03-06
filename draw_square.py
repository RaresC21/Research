import sys
#sys.path.append('/home/nbuser/library/')

from Polytope import Polytope
from Plane import Plane, Options
from symmetric import *
from far_point import *
from make_cube import *

import matplotlib.pyplot as plt
import numpy as np

#%matplotlib inline

plt.figure(figsize=(20,10))
plt.gca().set_aspect('equal', adjustable='box')

dim = 2
draw_circles = True

# cube = make_rotated_cube(dim)
# cube = make_cube(dim)
delta = 0.1
radius = np.sqrt(dim) * delta
depth = 0

options = Options(reflection = "sphere", \
                        delta = delta, radius = radius, depth = depth)
cube = make_cube(dim, options)
# zero = get_far_point(cube)
# print(np.linalg.norm(zero))
# cube.point = zero
# cube.point = np.zeros(dim)

centroid, lines = cube.approximate_centroid(iterations = 10000)
# lines = cube.cur_lines__

def plot_line(xs, vec, b):
    return b/vec[1] - (vec[0] * xs) / vec[1]

def plot_all_lines(polytope):
    for p in polytope.planes:
        xs = np.arange(-3,3, 0.1)
        ys = plot_line(xs, p.perp_vec, p.b)
        plt.plot(xs, ys)

def circ(x, y, r, phi):
  return (x + r * np.cos(phi), y + r*np.sin(phi))

def plot_circle(x, y, r):
    phis = np.arange(0,6.28,0.01)
    plt.plot(*circ(x, y, r, phis), c = 'r')

if draw_circles:
    for i in range(-1*int(1/delta), int(1/delta)+1):
        plot_circle(-1 - depth, i*delta, radius)
        plot_circle(1 + depth, i*delta, radius)
        plot_circle(i*delta, 1 + depth, radius)
        plot_circle(i*delta, -1 - depth, radius)

xx = 0
yy = 1
col = ['blue', 'red', 'black', 'green', 'purple']
c = 0
# for xx in range(0,dim):
#     for yy in range(xx+1,dim) :
pts = []
for i in range(len(lines)):
    pt1 = [lines[i][0][xx], lines[i][0][yy]]
    pt2 = [lines[i][1][xx], lines[i][1][yy]]
    # pts.append(lines[i][0])
    plt.plot([pt1[0], pt2[0]], [pt1[1], pt2[1]], color = col[c])
c += 1

plt.scatter(centroid[0], centroid[1], color = 'blue')

print("distance from centroid", np.linalg.norm(centroid))
print("centroid approximation", centroid)

# plt.xlim(-2,2)
# plt.ylim(-2,2)
# plt.xticks(np.arange(-1, 2, 1))
# plt.yticks(np.arange(-1, 2, 1))
plt.show()
