from Polytope import Polytope 

import numpy as np
from scipy.stats import ortho_group

def make_cube(dim, options, long = False):
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

    if long: 
        b[0] = 1000000
        b[1] = 1000000
        
    A = np.array(A)
    b = np.array(b)
    zero = np.array(zero)
    b = b * side / 2

    return Polytope(A, b, np.zeros(dim), options = options)    

def make_hyperplanes(dim, options) : 
    v1 = np.zeros(dim) 
    v2 = np.zeros(dim)
    v1[0] = 1
    v2[0] = -1
    b = np.array([1,1])
    A = np.array([v1, v2])
    return Polytope(A, b, np.zeros(dim), options = options)

def make_rotated_cube(dim, options): 
    A = ortho_group.rvs(dim = dim)
    A = np.vstack([A, -A])
    
    b = np.ones(2*dim)
    return Polytope(A, b, np.zeros(dim), options = options)
    