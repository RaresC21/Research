import matplotlib.pyplot as plt
import numpy as np

from Polytope import Polytope
from symmetric import *
from far_point import *
from make_cube import *

dim = 10
polytope = symmetric(dim)
# polytope = make_cube(dim)
polytope.point = get_far_point(polytope)
print("Initial Distance", np.linalg.norm(polytope.point))

see = [1000, 10000, 100000]

for s in range(len(see)):
    if s == 0: cc = False
    else: cc = True
    lines = polytope.reflection_lines(see[s], continue_ = cc)

    xs = [x for x in range(len(polytope.centroid_distance))]
    ys = polytope.centroid_distance
    dc = polytope.centroid_change

    # average_change = []
    # for i in range(len(dc)):
    #     if i < 100:
    #         average_change.append(sum(dc[:i]) / (i+1))
    #     else :
    #         average_change.append(sum(dc[i-100:i]) / 100)

    amnt = 0
    for i in range(len(ys) - 1):
        if ys[i] > ys[i+1]:
            amnt += 1

    print("answer", amnt, len(xs), amnt * 1.0 / len(xs))
    plt.plot(xs, ys)
    # plt.plot(xs, average_change)
    # plt.plot(xs, [1/(x+1) for x in xs])
    plt.show()
