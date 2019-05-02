import matplotlib.pyplot as plt
import numpy as np
import polytope as pc

from Polytope import Polytope
from Plane import Plane

from Utils import get_ratio

class Volume:
    def __init__(self, polytope):
        self.polytope = polytope

    def cuboid_volume(self, polytope):
        volume = 1
        length = [0] * polytope.dimension
        for plane in polytope.planes:
            nonzero = 0
            for indx, val in np.ndenumerate(plane.perp_vec):
                if abs(abs(val) - 1) <= 1e-8:
                    length[indx[0]] += plane.b
                    nonzero += 1
            if nonzero != 1: return -1
        for l in length:
            volume *= abs(l)
        return volume

    def compute(self):
        return self.approximate_volume(self.polytope, \
                                       axis = 0, \
                                       counter = 0, \
                                       looped = False)

    def approximate_volume(self, polytope, axis, counter, looped):
        cuboid = self.cuboid_volume(polytope)
        if cuboid > 0 and looped:
            return cuboid

        centroid, lines = polytope.approximate_centroid()

        cut = np.zeros(polytope.dimension)
        sgn = 1
        if counter < polytope.dimension: cut[axis] = 1
        else:
            cut[axis] = -1
            sgn = -1

        vol_ratio, next_pt = get_ratio(lines, Plane(cut, sgn * centroid[axis]))

        A = np.vstack([polytope.A, cut])
        Aprime = np.vstack([polytope.A, -1 * cut])
        b = np.append(polytope.b, sgn * centroid[axis])
        new_polytope = Polytope(A, b, centroid)
        new_polytope = new_polytope.remove_redundant()
        new_polytope.point = next_pt

        Aprime = np.vstack([polytope.A, -1 * cut])
        bprime = np.append(polytope.b, -1 * sgn * centroid[axis])
        other_half = Polytope(Aprime, bprime, centroid)

        exact = 0
        for i in range(10):
            exact += pc.volume(pc.Polytope(new_polytope.A, new_polytope.b)) / \
                     pc.volume(pc.Polytope(other_half.A, other_half.b))
        exact /= 10
        print("(real_ratio - approx_ratio)/ exact_ratio", (exact - vol_ratio) / exact)

        vol = self.approximate_volume(new_polytope,
                                          (axis + 1) % polytope.dimension,
                                          (counter + 1) % (2 * polytope.dimension),
                                          looped or ((counter + 1) == 2 * polytope.dimension))
        v = vol / vol_ratio
        return v + vol
