import numpy as np 
import networkx as nx
from statistics import mean, pstdev
from scipy.sparse import coo_matrix, csr_matrix, diags
from scipy.sparse.linalg import lsqr



def hierarchical_levels(graph, weight):
    """Returns the hierarchical levels of the nodes of a network in the ordering of graph.nodes().
    
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
    arXiv preprint arXiv:1908.04358.   
    """
    A = nx.adjacency_matrix(graph, weight=weight).transpose()
    k_in = csr_matrix(A.sum(axis=1))
    D_in = diags(A.sum(axis=1).A1, 0)
    L_in = D_in - A 
    return lsqr(L_in, k_in.toarray())[0]



def hierarchical_differences(graph, weight):
    """Returns the non-zero hierarchical differences over the edges of a network in the form of an adjacency matrix with the 
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
    arXiv preprint arXiv:1908.04358.   
    """
    A = nx.adjacency_matrix(graph, weight=weight).transpose()
    s = hierarchical_levels(graph, weight=weight)
    
    TD = []
    for i, j in zip(A.nonzero()[0], A.nonzero()[1]):
        TD.append(s[i].item() - s[j].item())
    TDmatrix = coo_matrix((TD, (A.nonzero()[0], A.nonzero()[1])))
    
    return TDmatrix.toarray(), mean(TD), pstdev(TD)



def hierarchical_coherence(graph, weight):
    """Returns the hierarchical levels, hierarchical differences over the edges of a network in the form of an adjacency matrix with the 
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
    arXiv preprint arXiv:1908.04358.   
    """
    A = nx.adjacency_matrix(graph,weight=weight).transpose()
    s = hierarchical_levels(graph, weight=weight)
    
    TD = []
    for i, j in zip(A.nonzero()[0], A.nonzero()[1]):
        TD.append(s[i].item() - s[j].item())
    TDmatrix = coo_matrix((TD, (A.nonzero()[0], A.nonzero()[1])))
    
    return s, TDmatrix.toarray(), mean(TD), pstdev(TD)


# returns a measure of equitable controllability over the full graph/network
def democracy_coefficient(graph, weight):
    """Returns the democracy coeffcient of a graph, a topological metric
    
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
    arXiv preprint arXiv:1908.04358.   
    """
    A = nx.adjacency_matrix(graph, weight=weight).transpose()
    s = hierarchical_levels(graph, weight=weight)
    TD = []
    for i, j in zip(A.nonzero()[0], A.nonzero()[1]):
        TD.append(s[i].item() - s[j].item())
    return 1 - mean(TD)



def influence_centrality(graph, weight, node):
    """Returns the influence centrality of a given node in the network with node labelling determined by graph.nodes().
    
    Parameters
    ----------
    graph : graph
       A NetworkX graph
       
    weight :  string or None
        If you have no weighted edges insert weight=None, otherwise weight='string', where string is your underlying weight attribute
    
    node : number or string
        Label of the node as determined by the makeup of the graph.nodes() call
        
    References
    ----------
    .. [1] Moutsinas, G., Shuaib, C., Guo, W., & Jarvis, S. (2019). 
    Graph hierarchy and spread of infections. 
    arXiv preprint arXiv:1908.04358.   
    """
    masked = np.ma.masked_equal(hierarchical_differences(graph, weight=weight)[0], 0)
    nodes = [x[0] for x in graph.nodes.items()]
    
    if node not in nodes:
        return print('This node is not part of the graph')
    try: 
        masked[node]
    except IndexError:
        return 1.0
    else:
        return 1.0 - np.ma.filled(masked[node].mean(), 0)


def total_influence_centrality(graph, weight):
    """Returns the influence centrality of all nodes in the network with node labelling determined by graph.nodes().
    
    Parameters
    ----------
    graph : graph
       A NetworkX graph
       
    weight :  string or None
        If you have no weighted edges insert weight=None, otherwise weight='string', where string is your underlying weight attribute
    
    node : number or string
        Label of the node as determined by the makeup of the graph.nodes() call
        
    References
    ----------
    .. [1] Moutsinas, G., Shuaib, C., Guo, W., & Jarvis, S. (2019). 
    Graph hierarchy and spread of infections. 
    arXiv preprint arXiv:1908.04358.   
    """
    nodes = [x[0] for x in graph.nodes.items()]
    ic = []
    for node in nodes:
        ic.append(influence_centrality(graph, weight, node))
    return ic
    
