import matplotlib.pyplot as plt
import numpy as np

from Polytope import Polytope
from symmetric import *
from far_point import *

def separation(polytope, iters): 
    lines = polytope.reflection_lines(iters)

    vec = np.random.random(polytope.dimension)
    vec /= np.linalg.norm(vec)
    
    total = np.zeros(polytope.dimension)
    f = []
    side = []
    for line in lines:
        pt = line[0]
        total += pt
        f.append(np.dot(vec, total) / (len(f) + 1))
        if np.dot(vec, pt) > 0: side.append(1)
        else: side.append(0)
        
        #print(len(f), f[-1])
        
    plt.plot([0,len(f)], [0,0])
    plt.plot(range(len(f)), f)
    plt.show()

    return side

dim = 50
poly = symmetric(dim)
zero = get_far_point(poly)
poly.point = zero

side = separation(poly, 10000)

a = []
b = []
cnt = 1
for i in range(len(side) - 1):
    if side[i] == 1 and side[i + 1] == 0:
        a.append(cnt)
        cnt = 0
    elif side[i] == 0 and side[i + 1] == 1:
        b.append(cnt)
        cnt = 0
    cnt += 1
print(a)
print(b)
print(np.mean(a))
print(np.mean(b))
