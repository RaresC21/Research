import matplotlib.pyplot as plt
import numpy as np

from Polytope import Polytope
from symmetric import *
from far_point import *
from make_cube import *

import collections

dim = 2
polytope = make_cube(dim) #symmetric(dim)

polytope.point = get_far_point(polytope)
polytope.reflection_lines(5000)

#print(np.mean(polytope.angles)) ## about 1/(sqrt(dim/2))

#plt.hist(x = polytope.angles, bins='auto')
#plt.show()

#facets = polytope.facet_hit
#plt.scatter(range(len(facets)), facets)
#plt.show()

def plot_facets_hit(data, a, b):
    bins = np.arange(min(data), max(data), 1)
    plt.xlim([min(data)-1, max(data)+1])

    d = collections.Counter(data[a:b]).values()
    #print(len(data), a, b, data[a:b], d, sum(data[a:b]))
    
    d = [x * 1.0 / sum(data[a:b]) for x in d]
    plt.bar([x for x in range(len(d))], d)
    # plt.show()

    return d

d1 = plot_facets_hit(polytope.facet_hit, 1000, 5000)
plt.show()
# d2 = plot_facets_hit(polytope.facet_hit, 1000, 10000)
# d3 = plot_facets_hit(polytope.facet_hit, 1000, 50000)
# d4 = plot_facets_hit(polytope.facet_hit, 1000, 1000000)
#
# print(d1)
# print(d2)
# print(d3)
# print(d4)
#
# diff = 0
# for a,b in zip(d1, d4):
#     diff += abs(a - b)
# print("5e3 vs 1e5", diff)
#
# diff = 0
# for a,b in zip(d2, d4):
#     diff += abs(a - b)
# print("1e4 vs 1e5", diff)
#
# diff = 0
# for a,b in zip(d3, d4):
#     diff += abs(a - b)
# print("5e4 vs 1e5", diff)
