
from Plane import Plane, Options 
from Polytope import Polytope

import numpy as np
import matplotlib.pyplot as plt

def symmetric(dim, num, options) :
    A = np.random.normal(0, 1, (num, dim))
    b = np.ones(num)

    for i in range(num//2, num):
        A[i] = -1 * A[i - (num//2)]
    zero = np.zeros(dim)

    ''' 
    angles = []
    for i in range(len(poly.planes)):
        for k in range(i + 1, len(poly.planes)):
            if abs(np.dot(poly.planes[i].perp_vec, poly.planes[k].perp_vec)) < 1 - 1e-8:
                angles.append(np.arccos(np.dot(poly.planes[i].perp_vec, poly.planes[k].perp_vec)))
    print(np.mean(angles))
    '''
    
    return Polytope(A, b, zero, options = options)