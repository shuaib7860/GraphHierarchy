from numpy import zeros, ones
from networkx import adjacency_matrix
from scipy.sparse import csr_matrix, diags, lil_matrix
from scipy.sparse.linalg import lsqr



def forward_hierarchical_levels(graph, weight):
    """Returns the forward hierarchical levels of the nodes of a network as an array.
    
    Parameters
    ----------
    graph : graph
       A NetworkX graph
       
    weight :  string or None
        If you have no weighted edges insert weight=None, otherwise weight='string', where string is your underlying weight attribute
    
    Returns
    -------
    forward hierarchical levels : array
        A Nx1 dimensional array indexed by the nodes, in the same order as graph.nodes, holding the value of their forward hierarchical levels.
    
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




def backward_hierarchical_levels(graph, weight):
    """Returns the backward hierarchical levels of the nodes of a network as an array.
    
    Parameters
    ----------
    graph : graph
       A NetworkX graph
       
    weight :  string or None
        If you have no weighted edges insert weight=None, otherwise weight='string', where string is your underlying weight attribute
    
    Returns
    -------
    backward hierarchical levels : array
        A Nx1 dimensional array indexed by the nodes, in the same order as graph.nodes, holding the value of their backward hierarchical levels.
    
    References
    ----------
    .. [1] Moutsinas, G., Shuaib, C., Guo, W., & Jarvis, S. (2019). 
    Graph hierarchy and spread of infections. 
    arXiv preprint arXiv:1908.04358."""
    
    A = adjacency_matrix(graph, weight=weight)
    k_in = csr_matrix(A.sum(axis=1))
    D_in = diags(A.sum(axis=1).A1, 0)
    L_in = D_in - A
    return lsqr(L_in, k_in.toarray())[0]



def hierarchical_levels(graph, weight):
    """Returns the hierarchical levels of the nodes of a network as an array which aids visualisation of the hierarchical structure in the network.
    
    Parameters
    ----------
    graph : graph
       A NetworkX graph
       
    weight :  string or None
        If you have no weighted edges insert weight=None, otherwise weight='string', where string is your underlying weight attribute
    
    Returns
    -------
    hierarchical levels : array
        A Nx1 dimensional array indexed by the nodes, in the same order as graph.nodes, holding the value of their hierarchical levels.
    
    References
    ----------
    .. [1] Moutsinas, G., Shuaib, C., Guo, W., & Jarvis, S. (2019). 
    Graph hierarchy and spread of infections. 
    arXiv preprint arXiv:1908.04358."""
    
    return forward_hierarchical_levels(graph, weight) - backward_hierarchical_levels(graph, weight)





def forward_hierarchical_differences(graph, weight):
    """Returns the forward hierarchical differences over the edges of a network in the form of a weighted adjacency matrix
    
    Parameters
    ----------
    graph : graph
       A NetworkX graph
       
    weight :  string or None
        If you have no weighted edges insert weight=None, otherwise weight='string', where string is your underlying weight attribute
    
    Returns
    -------
    forward hierarchical differences : array
        A NxN dimensional array representing a weighted adjancency matrix, with the edge weights corresponding to the forward hierarchical differences. 
        The column index represents the source node of the edge and the row index represents the destination node of the edge.
        
    References
    ----------
    .. [1] Moutsinas, G., Shuaib, C., Guo, W., & Jarvis, S. (2019). 
    Graph hierarchy and spread of infections. 
    arXiv preprint arXiv:1908.04358."""
    
    A = adjacency_matrix(graph, weight=weight).transpose()
    s = forward_hierarchical_levels(graph, weight=weight)
    TD = lil_matrix(A.shape, dtype=float)
    
    for i, j in zip(A.nonzero()[0], A.nonzero()[1]):
        TD[i,j] = s[i] - s[j]
   
    return TD.toarray()




def backward_hierarchical_differences(graph, weight):
    """Returns the backward hierarchical differences over the edges of a network in the form of a weighted adjacency matrix
    
    Parameters
    ----------
    graph : graph
       A NetworkX graph
       
    weight :  string or None
        If you have no weighted edges insert weight=None, otherwise weight='string', where string is your underlying weight attribute
    
    Returns
    -------
    backward hierarchical differences : array
        A NxN dimensional array representing a weighted adjancency matrix, with the edge weights corresponding to the backward hierarchical differences. 
        The column index represents the source node of the edge and the row index represents the destination node of the edge.
        
    References
    ----------
    .. [1] Moutsinas, G., Shuaib, C., Guo, W., & Jarvis, S. (2019). 
    Graph hierarchy and spread of infections. 
    arXiv preprint arXiv:1908.04358."""
    
    A = adjacency_matrix(graph, weight=weight)
    s = backward_hierarchical_levels(graph, weight=weight)
    TD = lil_matrix(A.shape, dtype=float)
    
    for i, j in zip(A.nonzero()[0], A.nonzero()[1]):
        TD[i,j] = s[i] - s[j]
   
    return TD.toarray()




def sparse_forward_hierarchical_differences(graph, weight):
    ''' Just a copy of the forward hierarchical differences function that returns the sparse matrix instead of the dense representation'''
    
    A = adjacency_matrix(graph, weight=weight).transpose()
    s = forward_hierarchical_levels(graph, weight=weight)
    TD = lil_matrix(A.shape, dtype=float)
    
    for i, j in zip(A.nonzero()[0], A.nonzero()[1]):
        TD[i,j] = s[i] - s[j]
   
    return TD



def sparse_backward_hierarchical_differences(graph, weight):
    ''' Just a copy of the backward hierarchical differences function that returns the sparse matrix instead of the dense representation'''
    
    A = adjacency_matrix(graph, weight=weight)
    s = backward_hierarchical_levels(graph, weight=weight)
    TD = lil_matrix(A.shape, dtype=float)
    
    for i, j in zip(A.nonzero()[0], A.nonzero()[1]):
        TD[i,j] = s[i] - s[j]
   
    return TD




def sparse_matrix_mean(sparse_matrix):
    '''A mean calculation for sparse matrices that does not count the zero elements. 
       So it calculates the mean of the hierarchical differences over the all incoming edges'''
    
    if sparse_matrix.sum() == 0:
        return 0
    else:
        return (sparse_matrix.sum()) / float(sparse_matrix.getnnz())




def forward_hierarchical_incoherence(graph, weight):
    """Returns the forward hierarchical differences over the edges of a network in the form of a weighted adjacency matrix,
    mean of the distribution of differences and standard deviation of this distribution.
    
    Parameters
    ----------
    graph : graph
       A NetworkX graph
       
    weight :  string or None
        If you have no weighted edges insert weight=None, otherwise weight='string', where string is your underlying weight attribute
    
    Returns
    -------
    forward hierarchical differences : array
        A NxN dimensional array representing a weighted adjancency matrix, with the edge weights corresponding to the forward hierarchical differences. 
        The column index represents the source node of the edge and the row index represents the destination node of the edge. 
        
    mean hierarchical difference : float
        The mean of the distribution of forward hierarchical differences.
    
    forward hierarchical incoherence : float
        The standard deviation of the distribution of forward hierarchical differences.
        
    References
    ----------
    .. [1] Moutsinas, G., Shuaib, C., Guo, W., & Jarvis, S. (2019). 
    Graph hierarchy and spread of infections. 
    arXiv preprint arXiv:1908.04358."""

    TD = sparse_forward_hierarchical_differences(graph, weight=weight)
    m = sparse_matrix_mean(TD.tocsr())
    std =  (sparse_matrix_mean(TD.power(2).tocsr()) - m**2)**0.5
    return TD.toarray(), m, std




def backward_hierarchical_incoherence(graph, weight):
    """Returns the backward  hierarchical differences over the edges of a network in the form of a weighted adjacency matrix,
    mean of the distribution of differences and standard deviation of this distribution.
    
    Parameters
    ----------
    graph : graph
       A NetworkX graph
       
    weight :  string or None
        If you have no weighted edges insert weight=None, otherwise weight='string', where string is your underlying weight attribute
    
    Returns
    -------
    backward hierarchical differences : array
        A NxN dimensional array representing a weighted adjancency matrix, with the edge weights corresponding to the backward hierarchical differences. 
        The column index represents the source node of the edge and the row index represents the destination node of the edge. 
        
    mean hierarchical difference : float
        The mean of the distribution of backward hierarchical differences.
    
    backward hierarchical incoherence : float
        The standard deviation of the distribution of backward hierarchical differences.
        
    References
    ----------
    .. [1] Moutsinas, G., Shuaib, C., Guo, W., & Jarvis, S. (2019). 
    Graph hierarchy and spread of infections. 
    arXiv preprint arXiv:1908.04358."""

    TD = sparse_backward_hierarchical_differences(graph, weight=weight)
    m = sparse_matrix_mean(TD.tocsr())
    std =  (sparse_matrix_mean(TD.power(2).tocsr()) - m**2)**0.5
    return TD.toarray(), m, std






# Returns a measure of equitable controllability over the full graph/network
def forward_democracy_coefficient(graph, weight):
    """Returns the forward democracy coeffcient of a graph, a topological network metric.
    
    Parameters
    ----------
    graph : graph
       A NetworkX graph
       
    weight :  string or None
        If you have no weighted edges insert weight=None, otherwise weight='string', where string is your underlying weight attribute
    
    Returns
    -------
    forward democracy coefficient : float
        forward democracy coefficient of a graph
        
    References
    ----------
    .. [1] Moutsinas, G., Shuaib, C., Guo, W., & Jarvis, S. (2019). 
    Graph hierarchy and spread of infections. 
    arXiv preprint arXiv:1908.04358."""
    
    TD = sparse_forward_hierarchical_differences(graph, weight=weight)
    return 1 - sparse_matrix_mean(TD.tocsr())




# Returns a measure of equitable controllability over the full graph/network
def backward_democracy_coefficient(graph, weight):
    """Returns the backward democracy coeffcient of a graph, a topological network metric.
    
    Parameters
    ----------
    graph : graph
       A NetworkX graph
       
    weight :  string or None
        If you have no weighted edges insert weight=None, otherwise weight='string', where string is your underlying weight attribute
    
    Returns
    -------
    backward democracy coefficient : float
        backward democracy coefficient of a graph
        
    References
    ----------
    .. [1] Moutsinas, G., Shuaib, C., Guo, W., & Jarvis, S. (2019). 
    Graph hierarchy and spread of infections. 
    arXiv preprint arXiv:1908.04358."""
    
    TD = sparse_backward_hierarchical_differences(graph, weight=weight)
    return 1 - sparse_matrix_mean(TD.tocsr())





def forward_influence_centrality(graph, weight):
    """Returns the forward influence centrality of the nodes in a network as an array.
    
    Parameters
    ----------
    graph : graph
       A NetworkX graph
       
    weight :  string or None
        If you have no weighted edges insert weight=None, otherwise weight='string', where string is your underlying weight attribute.
    
    Returns
    -------
    forward influence centrality : array
        A Nx1 dimensional array indexed by the nodes, in the same order as graph.nodes, holding the value of their forward influence centralities.
        
    References
    ----------
    .. [1] Moutsinas, G., Shuaib, C., Guo, W., & Jarvis, S. (2019). 
    Graph hierarchy and spread of infections. 
    arXiv preprint arXiv:1908.04358."""
    
    TD = sparse_forward_hierarchical_differences(graph, weight=weight)
    m = zeros((TD.shape[0], 1))
    for i in range(m.shape[0]):
        m[i] = sparse_matrix_mean(TD[i].tocsr())
    return ones((m.shape[0], 1)) - m



def backward_influence_centrality(graph, weight):
    """Returns the backward influence centrality of the nodes in a network as an array.
    
    Parameters
    ----------
    graph : graph
       A NetworkX graph
       
    weight :  string or None
        If you have no weighted edges insert weight=None, otherwise weight='string', where string is your underlying weight attribute.
    
    Returns
    -------
    backward influence centrality : array
        A Nx1 dimensional array indexed by the nodes, in the same order as graph.nodes, holding the value of their backward influence centralities.
        
    References
    ----------
    .. [1] Moutsinas, G., Shuaib, C., Guo, W., & Jarvis, S. (2019). 
    Graph hierarchy and spread of infections. 
    arXiv preprint arXiv:1908.04358."""
    
    TD = sparse_backward_hierarchical_differences(graph, weight=weight)
    m = zeros((TD.shape[0], 1))
    for i in range(m.shape[0]):
        m[i] = sparse_matrix_mean(TD[i].tocsr())
    return ones((m.shape[0], 1)) - m




def node_forward_influence_centrality(graph, weight, node):
    """Returns the forward influence centrality of the given node in the network.
    
    Parameters
    ----------
    graph : graph
       A NetworkX graph
       
    weight :  string or None
        If you have no weighted edges insert weight=None, otherwise weight='string', where string is your underlying weight attribute
    
    node : number
        Label of the node as determined by the indexing of the graph.nodes() call.
    
    Returns
    -------
    forward influence centrality : float
        A node's forward influence centrality.
        
    References
    ----------
    .. [1] Moutsinas, G., Shuaib, C., Guo, W., & Jarvis, S. (2019). 
    Graph hierarchy and spread of infections. 
    arXiv preprint arXiv:1908.04358."""
    
    return 1 - sparse_matrix_mean(sparse_forward_hierarchical_differences(graph, weight)[node].tocsr())




def node_backward_influence_centrality(graph, weight, node):
    """Returns the backward influence centrality of the given node in the network.
    
    Parameters
    ----------
    graph : graph
       A NetworkX graph
       
    weight :  string or None
        If you have no weighted edges insert weight=None, otherwise weight='string', where string is your underlying weight attribute
    
    node : number
        Label of the node as determined by the indexing of the graph.nodes() call.
    
    Returns
    -------
    backward influence centrality : float
        A node's backward influence centrality.
        
    References
    ----------
    .. [1] Moutsinas, G., Shuaib, C., Guo, W., & Jarvis, S. (2019). 
    Graph hierarchy and spread of infections. 
    arXiv preprint arXiv:1908.04358."""
    
    return 1 - sparse_matrix_mean(sparse_backward_hierarchical_differences(graph, weight)[node].tocsr())





def forward_hierarchical_metrics(graph, weight):
    ''' This function returns all the foundational node, edge and graph metrics a forward hierarchical/trophic approach yields.
    
    Parameters
    ----------
    graph : graph
       A NetworkX graph
       
    weight :  string or None
        If you have no weighted edges insert weight=None, otherwise weight='string', where string is your underlying weight attribute.
    
    Returns
    -------
    forward hierarchical levels : array
        A Nx1 dimensional array indexed by the nodes, in the same order as graph.nodes, holding the value of their forward hierarchical levels.
        
    forward influence centrality : array
        A Nx1 dimensional array indexed by the nodes, in the same order as graph.nodes, holding the value of their forward influence centralities.
    
    forward hierarchical differences : array
        A NxN dimensional array representing a weighted adjancency matrix, with the edge weights corresponding to the forward hierarchical differences. 
        The column index represents the source node of the edge and the row index represents the destination node of the edge.
     
    forward democracy coefficient : float
        forward democracy coefficient of a graph
    
    forward hierarchical incoherence : float
        The standard deviation of the distribution of forward hierarchical differences.
    
    References
    ----------
    .. [1] Moutsinas, G., Shuaib, C., Guo, W., & Jarvis, S. (2019). 
    Graph hierarchy and spread of infections. 
    arXiv preprint arXiv:1908.04358.
    '''
    
    A = adjacency_matrix(graph, weight=weight).transpose()
    s = forward_hierarchical_levels(graph, weight=weight)
    TD = lil_matrix(A.shape, dtype=float)
    
    for i, j in zip(A.nonzero()[0], A.nonzero()[1]):
        TD[i,j] = s[i] - s[j]
    
    m = sparse_matrix_mean(TD)  
    dc = 1 - m
    std = (sparse_matrix_mean(TD.power(2).tocsr()) - m**2)**0.5
    ic = forward_influence_centrality(graph, weight)
    
    return s, ic, TD, dc, std 





def backward_hierarchical_metrics(graph, weight):
    ''' This function returns all the foundational node, edge and graph metrics a backward hierarchical/trophic approach yields.
    
    Parameters
    ----------
    graph : graph
       A NetworkX graph
       
    weight :  string or None
        If you have no weighted edges insert weight=None, otherwise weight='string', where string is your underlying weight attribute.
    
    Returns
    -------
    backward hierarchical levels : array
        A Nx1 dimensional array indexed by the nodes, in the same order as graph.nodes, holding the value of their backward hierarchical levels.
        
    backward influence centrality : array
        A Nx1 dimensional array indexed by the nodes, in the same order as graph.nodes, holding the value of their backward influence centralities.
    
    backward hierarchical differences : array
        A NxN dimensional array representing a weighted adjancency matrix, with the edge weights corresponding to the backward hierarchical differences. 
        The column index represents the source node of the edge and the row index represents the destination node of the edge.
     
    backward democracy coefficient : float
        backward democracy coefficient of a graph
    
    backward hierarchical incoherence : float
        The standard deviation of the distribution of backward hierarchical differences.
    
    References
    ----------
    .. [1] Moutsinas, G., Shuaib, C., Guo, W., & Jarvis, S. (2019). 
    Graph hierarchy and spread of infections. 
    arXiv preprint arXiv:1908.04358.
    '''
    
    A = adjacency_matrix(graph, weight=weight)
    s = backward_hierarchical_levels(graph, weight=weight)
    TD = lil_matrix(A.shape, dtype=float)
    
    for i, j in zip(A.nonzero()[0], A.nonzero()[1]):
        TD[i,j] = s[i] - s[j]
    
    m = sparse_matrix_mean(TD)
    dc = 1 - m    
    std = (sparse_matrix_mean(TD.power(2).tocsr()) - m**2)**0.5
    ic = backward_influence_centrality(graph, weight)
    
    return s, ic, TD, dc, std 