import numpy as np
from Polytope import Polytope
from Plane import Plane, Options


def make_cube(dim, options):
    side = 2

    A = []
    b = []
    zero = [0] * dim

    for i in range(dim):
      vec = [0] * dim
      vec[i] = 1
      A.append(vec)

      vec1 = [0] * dim
      vec1[i] = -1
      A.append(vec1)

      b.append(1)
      b.append(1)

    A = np.array(A)
    b = np.array(b)
    zero = np.array(zero)
    b = b * side / 2

    return Polytope(A, b, zero, options = options)
