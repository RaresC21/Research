import numpy as np

def gram_schmidt(vectors):
    basis = []
    for v in vectors:
        w = v - np.sum( np.dot(v,b)*b  for b in basis )
        if (abs(w) > 1e-10).any():  
            basis.append(w/np.linalg.norm(w))
    return np.array(basis)


# given a vector, find an orthogonal extension forming a basis
def orthogonal_extend_basis(p): 
    dim = len(p)
    indx = -1
    
    for i, c in np.ndenumerate(p): 
        if abs(c) > 1e-8: 
            indx = i[0]
            break
        
    basis = [p]
    for i in range(dim): 
        if i == indx: 
            continue 
        cur = [0] * dim 
        cur[i] = 1
        basis.append(np.array(cur))
    return gram_schmidt(np.array(basis))