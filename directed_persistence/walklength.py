import numpy as np
from scipy.sparse.csgraph import shortest_path
from itertools import combinations, permutations
import gudhi


def walklength_filtration(digraph_data, max_dim=1):
    '''
    This function computes walklength filtration from a weighted directed graph.

    Args:
        digraph_data (numpy.array): Square 2d array where A[i,j] is the weight of directed edge (i,j); a zero value means there is no such edge.
        max_dim (int): maximum dimension to compute persistence; default is 1. Current version supports only max_dim <= 2.

    Returns:
        filt (gudhi.SimplexTree): simplicial filtration as a Gudhi simplex tree
    '''

    mat = shortest_path(digraph_data)

    vert_list = list(range(np.shape(mat)[0]))
    
    # Initialize filtration
    filt = gudhi.SimplexTree()

    ## Vertices (0-simplices)
    # Add all vertices with filtration value zero
    for v in vert_list:
        filt.insert([v,], 0)

    ## Edges (1-simplices)
    # Find filtration value for all edges
    for e in combinations(vert_list, 2):
        filt_val = min([mat[e[i-1],e[i]] for i in range(2)])
        filt.insert(e, filt_val)

    if max_dim > 0:
        ## Triangles (2-simplices)
        # Find filtration value for all triangles
        for t in combinations(vert_list, 3):
            clockwise = sum(sorted([mat[t[i-1],t[i]] for i in range(3)])[:-1])
            c_clockwise = sum(sorted([mat[t[i-2],t[i]] for i in range(3)])[:-1])
            # I'll explain later...
            filt_val = min(clockwise,c_clockwise)

            filt.insert(t, filt_val)

    if max_dim > 1:
        ## Tetrahedrons (3-simplices)
        # Find filtration value for all tetrahedrons
        d = 3
        for t in combinations(vert_list, d+1):
            possible_fvals = []
            for perm in permutations(range(d),d):
                cp = [*perm,d] # cyclic path
                min_path_cp = sum(sorted([mat[t[cp[i-1]],t[cp[i]]] for i in range(d+1)])[:-1])
                possible_fvals.append(min_path_cp)
            # I'll explain later...
            filt_val = min(possible_fvals)
            filt.insert(t, filt_val)

    if max_dim > 2:
        import warnings
        warnings.warn("Current version does not support homology of dimension higher than 2.")


    return filt
