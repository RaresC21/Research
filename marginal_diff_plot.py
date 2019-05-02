import matplotlib.pyplot as plt
import numpy as np

from Polytope import Polytope
from symmetric import *
from far_point import *

dim = 10
polytope = symmetric(dim)
polytope.point = get_far_point(polytope)
print("Initial Distance", np.linalg.norm(polytope.point))

centroid, lines = polytope.approximate_centroid()

#print(polytope.all_marginal_diffs)
#print(lines)
xs = [x * 50 for x in range(len(polytope.all_marginal_diffs))]
ys = polytope.all_marginal_diffs

plt.plot(xs, ys)
plt.show()

print("Final Distance", np.linalg.norm(centroid))
