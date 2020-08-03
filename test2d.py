import numpy as np

from Plane import Plane
from Polytope import Polytope, Options
from Volume import Volume

# arbitrary quadirlateral
# A = np.array([[2,1],[1,-1],[-3,4],[0,-1]])
# b = np.array([5,2,0,0])
A = np.array([[-1, 0], [0, -1],[1.0/100, 1]])
b = np.array([0,0,1])

delta = 0.01
radius = 0.1
depth = 0

options = Options(reflection = "facet", \
                        delta = delta, radius = radius, depth = depth)

triangle = Polytope(A, b, np.array([1.5,0.5]))

# area is still approximated, but basically exact
area = 50
# for i in range(250):
#     for k in range(150):
#         if triangle.pt_in_polytope(np.array([i,k]) / 100):
#             area += 1e-4

#area = side*side / 2
volume = Volume(triangle)
approx_area = volume.compute()

print("Approximate Area:", approx_area)
print("Exact Area:", area)
print("ratio:", approx_area / area)
