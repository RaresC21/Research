import matplotlib.pyplot as plt
import numpy as np
from Polytope import Polytope, Options

from centroid_change import * 

dim = 2
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
zero = np.array([0,-1])#np.array(zero)
b = b * side / 2

#right triangle
A = np.array([[-1,0],[1,1],[0,-1]])
b = np.array([0,2,2])
zero = np.array([0.5,0])

#sharp triangle
A = np.array([[-1,10],[50,1],[0,-1]])
b = np.array([2,2,2])
zero = np.array([-2,-1])

# hexagon
A = np.array([[1,1],[0,1],[0,-1],[1,-1],[-1,1],[-1,-1]])
b = np.array([4,4,4,4,8,8])
zero = np.array([0,0])

#right triangle
A = np.array([[-1,0],[1,1],[0,-1]])
b = np.array([0,2,2])
zero = np.array([0.5,0])

delta = 10
radius = np.sqrt(101)
depth = 10

cube = Polytope(A, b, zero, Options(reflection = "facet", \
                                    delta = delta, radius = radius, depth = depth))

draw_lines([[0,-2], [0,2], [4,-2], [0,-2]])

'''
plt.plot([-1,-1], [-1,1], color = 'black')
plt.plot([1,1], [-1,1], color = 'black')
plt.plot([-1,1], [1,1], color = 'black')
plt.plot([-1,1], [-1,-1], color = 'black')
'''

averages, centroids = cube.approximate_centroid(300, get_all = True)

xx = 0
yy = 1
i = 0.2
for j in range(len(averages[2:])):
  k = j + 2
  plt.plot([averages[k-1][xx], averages[k][xx]], \
           [averages[k-1][yy], averages[k][yy]], color = 'blue', alpha = i);
  i += .8 / len(centroids)

  plt.scatter(centroids[k][xx], centroids[k][yy], color = 'red', alpha = i)

#draw_lines([[0,-2], [0,2], [4,-2], [0,-2]])
  
'''
plt.plot([0,0], [-1,1], color = 'green')
plt.plot([-1,1], [0,0], color = 'green')
'''

plt.show()
