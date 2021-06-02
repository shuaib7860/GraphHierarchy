from numpy import zeros, ones, ndarray, average
from networkx import adjacency_matrix, Graph
from scipy.sparse import diags, lil_matrix, spmatrix
from scipy.sparse.linalg import lsqr



def forward_hierarchical_levels(graph, weight=None):
    """Returns the forward hierarchical levels of the nodes of a network as an array.
    
    Parameters
    ----------
    graph : Graph, array
       A NetworkX graph or numpy/sparse array
       
    weight :  string or None
        If you have weighted edges insert weight='string', where string is your underlying weight attribute. Only relevant if graph object is a networkx 
        graph instance. Otherwise the default is None.
    
    Returns
    -------
    forward hierarchical levels : array
        A Nx1 dimensional array indexed by the nodes, in the same order as graph.nodes if graph object, otherwise indexed in the same the numpy/sparse array, holding the value of their forward hierarchical levels.
    
    References
    ----------
    .. [1] Moutsinas, G., Shuaib, C., Guo, W., & Jarvis, S. (2019). 
    Graph hierarchy and spread of infections. 
    arXiv preprint arXiv:1908.04358."""
    
    if isinstance(graph, ndarray):
        A = graph.transpose()
        k_in = A.sum(axis=1)
        
    elif isinstance(graph, spmatrix):
        A = graph.transpose()
        k_in = A.sum(axis=1).A1
       
    elif isinstance(graph, Graph):
        A = adjacency_matrix(graph, weight=weight).transpose()
        k_in = A.sum(axis=1).A1
        
        
    D_in = diags(k_in, 0)
    L_in = D_in - A
    return lsqr(L_in, k_in)[0]





def backward_hierarchical_levels(graph, weight=None):
    """Returns the backward hierarchical levels of the nodes of a network as an array. This is the transpose of the original graph, so out-edges now become in-edges. 
    
    Parameters
    ----------
    graph : Graph, array
       A NetworkX graph or numpy/sparse array
       
    weight :  string or None
        If you have weighted edges insert weight='string', where string is your underlying weight attribute. Only relevant if graph object is a networkx 
        graph instance. Otherwise the default is None.
    
    Returns
    -------
    backward hierarchical levels : array
        A Nx1 dimensional array indexed by the nodes, in the same order as graph.nodes if graph object, otherwise indexed in the same the numpy/sparse array, holding the value of their forward hierarchical levels.
    
    References
    ----------
    .. [1] Moutsinas, G., Shuaib, C., Guo, W., & Jarvis, S. (2019). 
    Graph hierarchy and spread of infections. 
    arXiv preprint arXiv:1908.04358."""
    
    if isinstance(graph, ndarray):
        A = graph
        k_in = A.sum(axis=1)
        
    elif isinstance(graph, spmatrix):
        A = graph
        k_in = A.sum(axis=1).A1
       
    elif isinstance(graph, Graph):
        A = adjacency_matrix(graph, weight=weight)
        k_in = A.sum(axis=1).A1
        
        
    D_in = diags(k_in, 0)
    L_in = D_in - A
    return lsqr(L_in, k_in)[0]




def hierarchical_levels(graph, weight=None):
    """Returns the hierarchical levels of the nodes of a network as an array which aids visualisation of the hierarchical structure in the network.
    
    Parameters
    ----------
    graph : Graph
       A NetworkX graph or numpy/sparse array
       
    weight :  string or None
        If you have weighted edges insert weight='string', where string is your underlying weight attribute. Only relevant if graph object is a networkx 
        graph instance. Otherwise the default is None.
    
    Returns
    -------
    hierarchical levels : array
        A Nx1 dimensional array indexed by the nodes, in the same order as graph.nodes, holding the value of their hierarchical levels.
    
    References
    ----------
    .. [1] Moutsinas, G., Shuaib, C., Guo, W., & Jarvis, S. (2019). 
    Graph hierarchy and spread of infections. 
    arXiv preprint arXiv:1908.04358."""
    
    return 0.5*(forward_hierarchical_levels(graph, weight=weight) - backward_hierarchical_levels(graph, weight=weight))




def sparse_forward_hierarchical_differences(graph, weight=None):
    ''' Just a copy of the forward hierarchical differences function that returns the sparse matrix, instead of the dense representation, in lil format'''
    
    if isinstance(graph, (ndarray, spmatrix)):
        A = graph.transpose()
       
    elif isinstance(graph, Graph):
        A = adjacency_matrix(graph, weight=weight).transpose()
    
    s = forward_hierarchical_levels(graph, weight=weight)
    TD = lil_matrix(A.shape, dtype=float)
    
    for i, j in zip(A.nonzero()[0], A.nonzero()[1]):
        TD[i,j] = s[i] - s[j]
   
    return TD


def forward_hierarchical_differences(graph, weight=None):
    """Returns the forward hierarchical differences over the edges of a network in the form of a weighted adjacency matrix
    
    Parameters
    ----------
    graph : Graph, array
       A NetworkX graph or numpy/sparse array
       
    weight :  string
        If you have weighted edges insert weight='string', where string is your underlying weight attribute. Only relevant if graph object is a networkx 
        graph instance. Otherwise the default is None.
    
    Returns
    -------
    forward hierarchical differences : array
        A NxN dimensional array representing a weighted adjacency matrix, with the edge weights corresponding to the forward hierarchical differences. 
        The column index represents the source node of the edge and the row index represents the destination node of the edge.
        
    References
    ----------
    .. [1] Moutsinas, G., Shuaib, C., Guo, W., & Jarvis, S. (2019). 
    Graph hierarchy and spread of infections. 
    arXiv preprint arXiv:1908.04358."""
    
    TD = sparse_forward_hierarchical_differences(graph, weight=weight)
    return TD.toarray()





def sparse_backward_hierarchical_differences(graph, weight=None):
    ''' Just a copy of the backward hierarchical differences function that returns the sparse matrix, instead of the dense representation, in lil format'''
    
    if isinstance(graph, (ndarray, spmatrix)):
        A = graph
       
    elif isinstance(graph, Graph):
        A = adjacency_matrix(graph, weight=weight)
        
    s = backward_hierarchical_levels(graph, weight=weight)
    TD = lil_matrix(A.shape, dtype=float)
    
    for i, j in zip(A.nonzero()[0], A.nonzero()[1]):
        TD[i,j] = s[i] - s[j]
   
    return TD


def backward_hierarchical_differences(graph, weight=None):
    """Returns the backward hierarchical differences over the edges of a network in the form of a weighted adjacency matrix
    
    Parameters
    ----------
    graph : Graph, array
       A NetworkX graph or numpy/sparse array
       
    weight :  string or None
        If you have weighted edges insert weight='string', where string is your underlying weight attribute. Only relevant if graph object is a networkx 
        graph instance. Otherwise the default is None.
    
    Returns
    -------
    backward hierarchical differences : array
        A NxN dimensional array representing a weighted adjacency matrix, with the edge weights corresponding to the backward hierarchical differences. 
        The column index represents the source node of the edge and the row index represents the destination node of the edge.
        
    References
    ----------
    .. [1] Moutsinas, G., Shuaib, C., Guo, W., & Jarvis, S. (2019). 
    Graph hierarchy and spread of infections. 
    arXiv preprint arXiv:1908.04358."""
    
    TD = sparse_backward_hierarchical_differences(graph, weight=weight)
    return TD.toarray()





def forward_hierarchical_incoherence(graph, weight=None):
    """Returns the forward hierarchical differences over the edges of a network in the form of a weighted adjacency matrix,
    mean of the distribution of differences and standard deviation of this distribution.
    
    Parameters
    ----------
    graph : Graph, array
       A NetworkX graph or numpy/sparse array
       
    weight :  string or None
        If you have weighted edges insert weight='string', where string is your underlying weight attribute. Only relevant if graph object is a networkx 
        graph instance. Otherwise the default is None.
    
    Returns
    -------
    forward hierarchical differences : sparse array
        A NxN sparse dimensional sparse array representing a weighted adjancency matrix, with the edge weights corresponding to the forward hierarchical differences. 
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
    
    if isinstance(graph, ndarray):
        A = graph.transpose()
        TD = forward_hierarchical_differences(graph, weight=weight)
        m = average(TD, weights=A)
        m2 = average(TD**2, weights=A)
        
    elif isinstance(graph, spmatrix):
        A = graph.transpose()
        TD = sparse_forward_hierarchical_differences(graph, weight=weight).tocsc()
        m = (A.multiply(TD)).sum() / A.sum()
        m2 =  (A.multiply(TD.power(2))).sum() / A.sum()
    
    elif isinstance(graph, Graph):
        A = adjacency_matrix(graph, weight=weight).transpose()    
        TD = sparse_forward_hierarchical_differences(graph, weight=weight).tocsc()
        m = (A.multiply(TD)).sum() / A.sum()
        m2 =  (A.multiply(TD.power(2))).sum() / A.sum()
    
    std = (m2 - m**2)**0.5    
    return TD, m, std




def backward_hierarchical_incoherence(graph, weight=None):
    """Returns the backward  hierarchical differences over the edges of a network in the form of a weighted adjacency matrix,
    mean of the distribution of differences and standard deviation of this distribution.
    
    Parameters
    ----------
    graph : Graph, array
       A NetworkX graph or numpy/sparse array
       
    weight :  string
        If you have weighted edges insert weight='string', where string is your underlying weight attribute. Only relevant if graph object is a networkx 
        graph instance. Otherwise the default is None.
    
    Returns
    -------
    backward hierarchical differences : sparse array
        A NxN dimensional sparse array representing a weighted adjancency matrix, with the edge weights corresponding to the backward hierarchical differences. 
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

    if isinstance(graph, ndarray):
        A = graph
        TD = backward_hierarchical_differences(graph, weight=weight)
        m = average(TD, weights=A)
        m2 = average(TD**2, weights=A)
        
    elif isinstance(graph, spmatrix):
        A = graph
        TD = sparse_backward_hierarchical_differences(graph, weight=weight).tocsr()
        m = (A.multiply(TD)).sum() / A.sum()
        m2 =  (A.multiply(TD.power(2))).sum() / A.sum()
    
    elif isinstance(graph, Graph):
        A = adjacency_matrix(graph, weight=weight)  
        TD = sparse_backward_hierarchical_differences(graph, weight=weight).tocsr()
        m = (A.multiply(TD)).sum() / A.sum()
        m2 =  (A.multiply(TD.power(2))).sum() / A.sum()
    
    std = (m2 - m**2)**0.5    
    return TD, m, std





# Returns a measure of equitable controllability over the full graph/network
def forward_democracy_coefficient(graph, weight=None):
    """Returns the forward democracy coeffcient of a graph, a topological network metric.
    
    Parameters
    ----------
    graph : Graph, array
       A NetworkX graph or numpy/sparse array
       
    weight :  string
        If you have weighted edges insert weight='string', where string is your underlying weight attribute. Only relevant if graph object is a networkx 
        graph instance. Otherwise the default is None.
        
    Returns
    -------
    forward democracy coefficient : float
        forward democracy coefficient of a graph
        
    References
    ----------
    .. [1] Moutsinas, G., Shuaib, C., Guo, W., & Jarvis, S. (2019). 
    Graph hierarchy and spread of infections. 
    arXiv preprint arXiv:1908.04358."""
    
    if isinstance(graph, ndarray):
        A = graph.transpose()
        TD = forward_hierarchical_differences(graph, weight=weight)
        m = average(TD, weights=A)
        
    elif isinstance(graph, spmatrix):
        A = graph.transpose()
        TD = sparse_forward_hierarchical_differences(graph, weight=weight).tocsc()
        m = (A.multiply(TD)).sum() / A.sum()
        
    elif isinstance(graph, Graph):
        A = adjacency_matrix(graph, weight=weight).transpose()
        TD = sparse_forward_hierarchical_differences(graph, weight=weight).tocsc()
        m = (A.multiply(TD)).sum() / A.sum()
        
    return 1 - m




# Returns a measure of equitable controllability over the full graph/network
def backward_democracy_coefficient(graph, weight=None):
    """Returns the backward democracy coeffcient of a graph, a topological network metric.
    
    Parameters
    ----------
    graph : Graph, array
       A NetworkX graph or numpy/sparse array
       
    weight :  string or None
        If you have weighted edges insert weight='string', where string is your underlying weight attribute. Only relevant if graph object is a networkx 
        graph instance. Otherwise the default is None.
    
    Returns
    -------
    backward democracy coefficient : float
        backward democracy coefficient of a graph
        
    References
    ----------
    .. [1] Moutsinas, G., Shuaib, C., Guo, W., & Jarvis, S. (2019). 
    Graph hierarchy and spread of infections. 
    arXiv preprint arXiv:1908.04358."""
    
    if isinstance(graph, ndarray):
        A = graph
        TD = backward_hierarchical_differences(graph, weight=weight)
        m = average(TD, weights=A)
        
    elif isinstance(graph, spmatrix):
        A = graph
        TD = sparse_backward_hierarchical_differences(graph, weight=weight).tocsr()
        m = (A.multiply(TD)).sum() / A.sum()
        
    elif isinstance(graph, Graph):
        A = adjacency_matrix(graph, weight=weight)
        TD = sparse_backward_hierarchical_differences(graph, weight=weight).tocsr()
        m = (A.multiply(TD)).sum() / A.sum()
    
    return 1 - m



def node_forward_influence_centrality(graph, node, weight=None):
    """Returns the forward influence centrality of the given node in the network.
    
    Parameters
    ----------
    graph : Graph, array
       A NetworkX graph or numpy/sparse array
    
    node : number
        Label of the node as determined by the indexing of the graph.nodes() call or the index of the numpy/sparse array.
        
    weight :  string
        If you have weighted edges insert weight='string', where string is your underlying weight attribute. Only relevant if graph object is a networkx 
        graph instance. Otherwise the default is None.
        
    Returns
    -------
    forward influence centrality : float
        A node's forward influence centrality.
        
    References
    ----------
    .. [1] Moutsinas, G., Shuaib, C., Guo, W., & Jarvis, S. (2019). 
    Graph hierarchy and spread of infections. 
    arXiv preprint arXiv:1908.04358."""
    
    if isinstance(graph, ndarray):
        A = graph.transpose()
        index = node
        TD = forward_hierarchical_differences(graph, weight=weight) 
        if A[index].sum() == 0:
            m = 0
        else:
            m = average(TD[index], weights=A[index])
        
    elif isinstance(graph, spmatrix):
        A = graph.transpose()
        index = node
        TD = sparse_forward_hierarchical_differences(graph, weight=weight).tocsc()
        if A[index].sum() == 0:
            m = 0
        else:
            m = (A[index].multiply(TD[index])).sum() / A[index].sum()
            
    elif isinstance(graph, Graph):
        A = adjacency_matrix(graph, weight=weight).transpose()
        index = list(graph.nodes).index(node)
        TD = sparse_forward_hierarchical_differences(graph, weight=weight).tocsc()
        if A[index].sum() == 0:
            m = 0
        else:
            m = (A[index].multiply(TD[index])).sum() / A[index].sum()
    
    return 1 - m




def node_backward_influence_centrality(graph, node, weight=None):
    """Returns the backward influence centrality of the given node in the network.
    
    Parameters
    ----------
    graph : Graph array
       A NetworkX graph or numpy/sparse array
        
    node : number
        Label of the node as determined by the indexing of the graph.nodes() call or the index of the numpy/sparse array.

    weight :  string or None
        If you have weighted edges insert weight='string', where string is your underlying weight attribute. Only relevant if graph object is a networkx 
        graph instance, otherwise the default is None.
    
    
    Returns
    -------
    backward influence centrality : float
        A node's backward influence centrality.
        
    References
    ----------
    .. [1] Moutsinas, G., Shuaib, C., Guo, W., & Jarvis, S. (2019). 
    Graph hierarchy and spread of infections. 
    arXiv preprint arXiv:1908.04358."""
    
    if isinstance(graph, ndarray):
        A = graph
        index = node
        TD = backward_hierarchical_differences(graph, weight=weight) 
        if A[index].sum() == 0:
            m = 0
        else:
            m = average(TD[index], weights=A[index])
        
    elif isinstance(graph, spmatrix):
        A = graph
        index = node
        TD = sparse_backward_hierarchical_differences(graph, weight=weight).tocsr()
        if A[index].sum() == 0:
            m = 0
        else:
            m = (A[index].multiply(TD[index])).sum() / A[index].sum()
            
    elif isinstance(graph, Graph):
        A = adjacency_matrix(graph, weight=weight)
        index = list(graph.nodes).index(node)
        TD = sparse_backward_hierarchical_differences(graph, weight=weight).tocsr()
        if A[index].sum() == 0:
            m = 0
        else:
            m = (A[index].multiply(TD[index])).sum() / A[index].sum()
    
    return 1 - m




def forward_influence_centrality(graph, weight=None):
    """Returns the forward influence centrality of the nodes in a network as an array.
    
    Parameters
    ----------
    graph : Graph, array
       A NetworkX graph or numpy/sparse array
       
    weight :  string or None
        If you have weighted edges insert weight='string', where string is your underlying weight attribute. Only relevant if graph object is a networkx 
        graph instance, otherwise the default is None.    
    Returns
    -------
    forward influence centrality : array
        A Nx1 dimensional array indexed by the nodes, in the same order as graph.nodes, holding the value of their forward influence centralities.
        
    References
    ----------
    .. [1] Moutsinas, G., Shuaib, C., Guo, W., & Jarvis, S. (2019). 
    Graph hierarchy and spread of infections. 
    arXiv preprint arXiv:1908.04358."""
    
    if isinstance(graph, ndarray):
        A = graph.transpose()
        TD = forward_hierarchical_differences(graph, weight=weight)
        m = zeros((TD.shape[0], 1))
    
        for i in range(m.shape[0]):
            if A[i].sum() == 0:
                m[i] = 0
            else:
                m[i] = average(TD[i], weights=A[i])
        
    elif isinstance(graph, spmatrix):
        A = graph.transpose()
        TD = sparse_forward_hierarchical_differences(graph, weight=weight).tocsc()
        m = zeros((TD.shape[0], 1))
    
        for i in range(m.shape[0]):
            if A[i].sum() == 0:
                m[i] = 0
            else:
                m[i] = (A[i].multiply(TD[i])).sum() / A[i].sum()
        
    elif isinstance(graph, Graph):
        A = adjacency_matrix(graph, weight=weight).transpose()
        TD = sparse_forward_hierarchical_differences(graph, weight=weight).tocsc()
        m = zeros((TD.shape[0], 1))
    
        for i in range(m.shape[0]):
            if A[i].sum() == 0:
                m[i] = 0
            else:
                m[i] = (A[i].multiply(TD[i])).sum() / A[i].sum()
    
    return ones((m.shape[0], 1)) - m





def backward_influence_centrality(graph, weight=None):
    """Returns the backward influence centrality of the nodes in a network as an array.
    
    Parameters
    ----------
    graph : Graph, array
        A NetworkX graph or numpy/sparse array
       
    weight :  string
        If you have weighted edges insert weight='string', where string is your underlying weight attribute. Only relevant if graph object is a networkx 
        graph instance, otherwise the default is None.
        
    Returns
    -------
    backward influence centrality : array
        A Nx1 dimensional array indexed by the nodes, in the same order as graph.nodes, holding the value of their backward influence centralities.
        
    References
    ----------
    .. [1] Moutsinas, G., Shuaib, C., Guo, W., & Jarvis, S. (2019). 
    Graph hierarchy and spread of infections. 
    arXiv preprint arXiv:1908.04358."""
            
    if isinstance(graph, ndarray):
        A = graph
        TD = backward_hierarchical_differences(graph, weight=weight)
        m = zeros((TD.shape[0], 1))
    
        for i in range(m.shape[0]):
            if A[i].sum() == 0:
                m[i] = 0
            else:
                m[i] = average(TD[i], weights=A[i])
        
    elif isinstance(graph, spmatrix):
        A = graph
        TD = sparse_backward_hierarchical_differences(graph, weight=weight).tocsr()
        m = zeros((TD.shape[0], 1))
    
        for i in range(m.shape[0]):
            if A[i].sum() == 0:
                m[i] = 0
            else:
                m[i] = (A[i].multiply(TD[i])).sum() / A[i].sum()
        
    elif isinstance(graph, Graph):
        A = adjacency_matrix(graph, weight=weight)
        TD = sparse_backward_hierarchical_differences(graph, weight=weight).tocsr()
        m = zeros((TD.shape[0], 1))
    
        for i in range(m.shape[0]):
            if A[i].sum() == 0:
                m[i] = 0
            else:
                m[i] = (A[i].multiply(TD[i])).sum() / A[i].sum()
                
    return ones((m.shape[0], 1)) - m






def forward_hierarchical_metrics(graph, weight=None):
    ''' This function returns all the foundational node, edge and graph metrics a forward hierarchical approach yields.
    
    Parameters
    ----------
    graph : Graph, array
       A NetworkX graph or numpy/sparse array
       
    weight :  string or None
        If you have weighted edges insert weight='string', where string is your underlying weight attribute. Only relevant if graph object is a networkx 
        graph instance. Otherwise the default is None.
        
    Returns
    -------
    forward hierarchical levels : array
        A Nx1 dimensional array indexed by the nodes, in the same order as graph.nodes, holding the value of their forward hierarchical levels.
        
    forward influence centrality : array
        A Nx1 dimensional array indexed by the nodes, in the same order as graph.nodes, holding the value of their forward influence centralities.
    
    forward hierarchical differences : sparse array
        A NxN dimensional sparse array representing a weighted adjancency matrix, with the edge weights corresponding to the forward hierarchical differences. 
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
    s = forward_hierarchical_levels(graph, weight=weight)
    ic = forward_influence_centrality(graph, weight=weight)
    a = forward_hierarchical_incoherence(graph, weight=weight)
    
    return s, ic, a[0], 1 - a[1], a[2] 





def backward_hierarchical_metrics(graph, weight=None):
    ''' This function returns all the foundational node, edge and graph metrics a backward hierarchical approach yields.
    
    Parameters
    ----------
    graph : Graph, array
       A NetworkX graph or numpy/sparse array
       
    weight :  string or None
        If you have weighted edges insert weight='string', where string is your underlying weight attribute. Only relevant if graph object is a networkx 
        graph instance. Otherwise the default is None.
        
    Returns
    -------
    backward hierarchical levels : array
        A Nx1 dimensional array indexed by the nodes, in the same order as graph.nodes, holding the value of their backward hierarchical levels.
        
    backward influence centrality : array
        A Nx1 dimensional array indexed by the nodes, in the same order as graph.nodes, holding the value of their backward influence centralities.
    
    backward hierarchical differences : sparse array
        A NxN dimensional sparse array representing a weighted adjancency matrix, with the edge weights corresponding to the backward hierarchical differences. 
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
    
    s = backward_hierarchical_levels(graph, weight=weight)
    ic = backward_influence_centrality(graph, weight=weight)
    a = backward_hierarchical_incoherence(graph, weight=weight)
    
    return s, ic, a[0], 1 - a[1], a[2] 




