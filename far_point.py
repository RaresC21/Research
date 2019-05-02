import numpy as np
from Polytope import Polytope
from symmetric import *

def get_far_point(polytope):
    lines = polytope.reflection_lines(100)

    dist = 0
    pt = -1
    for l in lines:
        a = l[0]
        b = l[1]
        p1 = 0.9 * a + 0.1 * b
        d1 = np.linalg.norm(p1)

        p2 = 0.1 * a + 0.9 * b
        d2 = np.linalg.norm(p2)

        if d1 > dist:
            dist = d1
            pt = p1
        if d2 > dist:
            dist = d2
            pt = p2
    return pt
