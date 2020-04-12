from numpy import zeros, ones
from networkx import adjacency_matrix
from scipy.sparse import csr_matrix, diags
from scipy.sparse.linalg import lsqr


def hierarchical_levels(graph, weight):
    """Returns the hierarchical levels of the nodes of a network as an array in the ordering/indexing of graph.nodes().
    
    Parameters
    ----------
    graph : graph
       A NetworkX graph
       
    weight :  string or None
        If you have no weighted edges insert weight=None, otherwise weight='string', where string is your underlying weight attribute
           
    References
    ----------
    .. [1] Moutsinas, G., Shuaib, C., Guo, W., & Jarvis, S. (2019). 
    Graph hierarchy and spread of infections. 
    arXiv preprint arXiv:1908.04358."""
    A = adjacency_matrix(graph, weight=weight).transpose()
    k_in = csr_matrix(A.sum(axis=1))
    D_in = diags(A.sum(axis=1).A1, 0)
    L_in = D_in - A
    return lsqr(L_in, k_in.toarray())[0]


def hierarchical_differences(graph, weight):
    """Returns the non-zero hierarchical differences over the edges of a network in the form of an adjacency matrix with the 
    ordering/indexing of graph.nodes(). 
    
    Parameters
    ----------
    graph : graph
       A NetworkX graph
       
    weight :  string or None
        If you have no weighted edges insert weight=None, otherwise weight='string', where string is your underlying weight attribute
        
        
    References
    ----------
    .. [1] Moutsinas, G., Shuaib, C., Guo, W., & Jarvis, S. (2019). 
    Graph hierarchy and spread of infections. 
    arXiv preprint arXiv:1908.04358."""
    A = adjacency_matrix(graph, weight=weight).transpose()
    s = hierarchical_levels(graph, weight=weight)
    TD = A.copy()
    TD = TD.astype(float)
    
    for i, j in zip(TD.nonzero()[0], TD.nonzero()[1]):
        TD[i,j] = s[i] - s[j]
   
    return TD.toarray()


def sparse_hierarchical_differences(graph, weight):
    ''' Just a copy of the hierarchical differences function that returns the sparse matrix instead of the dense representation'''
    A = adjacency_matrix(graph, weight=weight).transpose()
    s = hierarchical_levels(graph, weight=weight)
    TD = A.copy()
    TD = TD.astype(float)
    
    for i, j in zip(TD.nonzero()[0], TD.nonzero()[1]):
        TD[i,j] = s[i] - s[j]
   
    return TD


def sparse_matrix_mean(sparse_matrix):
    '''A mean calculation for sparse matrices that does not count the zero elements'''
    if sparse_matrix.sum(dtype=float) == 0:
        return 0
    else:
        return (sparse_matrix.sum(dtype=float)) / float(sparse_matrix.count_nonzero())


def hierarchical_coherence(graph, weight):
    """Returns the hierarchical differences over the edges of a network in the form of an adjacency matrix with the 
    ordering/indexing of graph.nodes(), mean of the distribution of differences and standard deviation of this distribution.
    
    Parameters
    ----------
    graph : graph
       A NetworkX graph
       
    weight :  string or None
        If you have no weighted edges insert weight=None, otherwise weight='string', where string is your underlying weight attribute
        
        
    References
    ----------
    .. [1] Moutsinas, G., Shuaib, C., Guo, W., & Jarvis, S. (2019). 
    Graph hierarchy and spread of infections. 
    arXiv preprint arXiv:1908.04358."""
    TD = sparse_hierarchical_differences(graph, weight=weight)
    std =  (sparse_matrix_mean(TD.power(2)) - sparse_matrix_mean(TD)**2)**0.5
    return TD.toarray(), sparse_matrix_mean(TD), std


# Returns a measure of equitable controllability over the full graph/network
def democracy_coefficient(graph, weight):
    """Returns the democracy coeffcient of a graph, a topological metric.
    
    Parameters
    ----------
    graph : graph
       A NetworkX graph
       
    weight :  string or None
        If you have no weighted edges insert weight=None, otherwise weight='string', where string is your underlying weight attribute
        
        
    References
    ----------
    .. [1] Moutsinas, G., Shuaib, C., Guo, W., & Jarvis, S. (2019). 
    Graph hierarchy and spread of infections. 
    arXiv preprint arXiv:1908.04358."""
    TD = sparse_hierarchical_differences(graph, weight=weight)
    return 1 - sparse_matrix_mean(TD)


def influence_centrality(graph, weight):
    """Returns the influence centrality of all nodes in the network as an array in the ordering/indexing of graph.nodes().
    
    Parameters
    ----------
    graph : graph
       A NetworkX graph
       
    weight :  string or None
        If you have no weighted edges insert weight=None, otherwise weight='string', where string is your underlying weight attribute.
        
    References
    ----------
    .. [1] Moutsinas, G., Shuaib, C., Guo, W., & Jarvis, S. (2019). 
    Graph hierarchy and spread of infections. 
    arXiv preprint arXiv:1908.04358."""
    TD = sparse_hierarchical_differences(graph, weight=weight)
    m = zeros((TD.shape[0], 1))
    for i in range(TD.shape[0]):
        m[i] = sparse_matrix_mean(TD[i])
    return ones((TD.shape[0], 1)) - m


def node_influence_centrality(graph, weight, node):
    """Returns the influence centrality of a given node in the network with node indexing determined by graph.nodes().
    
    Parameters
    ----------
    graph : graph
       A NetworkX graph
       
    weight :  string or None
        If you have no weighted edges insert weight=None, otherwise weight='string', where string is your underlying weight attribute
    
    node : number or string
        Label of the node as determined by the makeup of the graph.nodes() call.
        
    References
    ----------
    .. [1] Moutsinas, G., Shuaib, C., Guo, W., & Jarvis, S. (2019). 
    Graph hierarchy and spread of infections. 
    arXiv preprint arXiv:1908.04358."""
    return 1 - sparse_matrix_mean(sparse_hierarchical_differences(graph, weight)[node])