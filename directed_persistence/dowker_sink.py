import numpy as np
from itertools import chain, combinations
import gudhi


def dowker_sink(A, max_dim=1):
    '''
    This function computes Dowker sink filtration from a network, that is, a complete weighted directed graph.

    Args:
        A (numpy.array): Square 2d array where A[i,j] is the weight from node i to node j; can be symmetric or asymmetric.
        max_dim (int): maximum dimension to compute persistence; default is 1.

    Returns:
        filt (gudhi.SimplexTree): simplicial filtration as a Gudhi simplex tree
    '''
    
    # Get list of all vertices
    vert_list = list(range(np.shape(A)[0]))

    # Initialize filtration
    filt = gudhi.SimplexTree()

    for v in vert_list:
        # add the vertex with filtration value equal to the min weight
        # of all edges going out of v
        filt.insert([v,], min(A[v,:]))

    # Loop over all possible simplices to find each filtration value
    for simp in list(chain.from_iterable(combinations(vert_list, r) for r in range(2,max_dim+3))):
        filt_value = np.min(np.amax([A[i,:] for i in simp], axis=0))
        filt.insert(simp, filt_value)


    return filt
