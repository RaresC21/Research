from collections import defaultdict
from scipy.optimize import linprog
from symmetric import *
from make_cube import *

def get_neighbors(polytope):
    neighbors = defaultdict(list)
    for i in range(polytope.num_planes):
        for k in range(i + 1, polytope.num_planes):
            A = polytope.A.copy()
            b = polytope.b.copy()

            A = np.vstack([A, -1 * A[i]])
            A = np.vstack([A, -1 * A[k]])
            
            b = np.append(b, -b[i])
            b = np.append(b, -b[k])

            res = linprog(np.ones(polytope.dimension), A_ub=A, b_ub=b, bounds=(None,None))
            if res.status == 0:
                neighbors[i].append(k)
                neighbors[k].append(i)
                
    return neighbors


p = symmetric(10)

adj = get_neighbors(p)
print(len(adj))
for a in adj:
    print(a, len(adj[a]))
