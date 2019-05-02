import numpy as np

def read_graph(): 

    # input graph in file "graph"
    inp = open("graph")
    V = int(inp.readline())
    E = int(inp.readline())
    
    out = [[] for i in range(V)]
    inn = [[] for i in range(V)]
    
    for i in range(E):
        edge = inp.readline().split()
        a = int(edge[0]) - 1
        b = int(edge[1]) - 1
        
        out[a].append(i)
        inn[b].append(i)
        
    rows = V + 2*E - 2
    cols = E 
    A = [[0 for x in range(cols)] for y in range(rows)] 
    b = [0 for x in range(rows)]
    
    for i in range(0,V-2):
        for e in inn[i+1]:
            A[i][e] = 1
        for e in out[i+1]:
            A[i][e] = -1
        
    e = 0
    for i in range(V-2, V+E-2):
        A[i][e] = 1
        e += 1

        b[i] = 1

    e = 0
    for i in range(V+E-2, V+2*E-2):
        A[i][e] = -1
        e += 1

    return np.array(A), np.array(b)

    
