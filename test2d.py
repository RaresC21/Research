import numpy as np

from Plane import Plane
from Polytope import Polytope
import Volume

# arbitrary quadirlateral
A = np.array([[2,1],[1,-1],[-3,4],[0,-1]])
b = np.array([5,2,0,0])

triangle = Polytope(A, b, np.array([1.5,0.5]))

# area is still approximated, but basically exact
area = 0
for i in range(250):
    for k in range(150):
        if triangle.pt_in_polytope(np.array([i,k]) / 100):
            area += 1e-4

#area = side*side / 2
approx_area = Volume.approximate_volume(triangle)

print("Approximate Area:", approx_area)
print("Exact Area:", area)
print("ratio:", approx_area / area)
