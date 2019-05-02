import numpy as np
import polytope as poly
import matplotlib.pyplot as plt

from Plane import Plane, Options
from Polytope import Polytope
from Volume import Volume

from flow_polytope import read_graph

def flow_test(dim, iterations, accuracy, options):
    A,b = read_graph()
    print("A = ", A)
    print("b = ", b)
    
    pp = poly.Polytope(A,b)
    c = poly.cheby_ball(pp)
    zero = c[1]
    
    polytope = Polytope(A, b, zero)
    volume = Volume(polytope).compute()

    exact = 0
    all_ = []
    for i in range(100):
        all_.append(poly.volume(poly.Polytope(A, b)))
        exact += all_[-1]
    exact /= 100

    print("exact volume:", exact)
    print("std. dev (of \"exact\" volume)", np.std(all_));
    print("my approx:", volume)
    return volume / exact

def test_cube(dim, iterations, accuracy, options) :
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
    
    cube = Polytope(A, b, zero, options)
    
    volume, it = Volume(cube).compute(iterations, accuracy)

    exact_volume = poly.volume(poly.Polytope(cube.A, cube.b))
    print("exact volume:", exact_volume)
    return volume / exact_volume, it


'''
options = Options(reflection = "facet", \
                  delta = 0.1, radius = 0.1 * np.sqrt(10), depth = 0)
dim = [2,3,5, 10, 15, 20]
iterations = []
for i in dim:
    ratio, it = test(dim = i, iterations = 10, accuracy = 0.1, options = options)
    iterations.append(np.mean(it))
print(dim, iterations)
plt.plot(dim, iterations)
plt.xlabel('dimension')
plt.ylabel('iterations')
plt.show()
'''

options = Options(reflection = "sphere", \
                  delta = 0.01, radius = 0.01 * np.sqrt(5), depth = 0)

ratio = flow_test(dim = 20, iterations = 1, accuracy = 0.05, options = options);

print("ratio (approx)/(real):", ratio)


