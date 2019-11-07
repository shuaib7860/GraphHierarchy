import numpy as np 
import networkx as nx
from statistics import mean, pstdev
from scipy.sparse import coo_matrix, csr_matrix, diags, 
from scipe.sparse.linalg import lsqr


#The index postion of the element in the vector matches the networkx label of the node, so 0th element of vector is 0th labelled node. 
def HierarchicalLevel(graph):
    A = nx.adjacency_matrix(graph).transpose()
    k_in = A.sum(axis=1)
    D_in = np.diagflat(k_in) #check to make sure this is a sparse matrix, it's not sparse
    L_in = D_in - A
    Lt_in = np.linalg.pinv(L_in)
    s = Lt_in * k_in
    return s
# maybe return s as a dictionary with the key as the label of the node and the value is the gtl
# Could implement linear solver for normal numpy array in above code


# This is the more efficient piece of code, I think, need to confirm; confirmed 
# k needs to be an array and so the csr function does have a toarray method but the normal numpy matrix object does not
def HierarchicalLevelSparse(graph):
    A = nx.adjacency_matrix(graph).transpose()
    k_in = csr_matrix(A.sum(axis=1))
    D_in = diags(A.sum(axis=1).A1, 0)
    L_in = D_in - A
    s = lsqr(L_in, k_in.toarray())
    return s

def HierarchicalDifferences(graph):
    A = nx.adjacency_matrix(graph).transpose()
    s = HierarchicalLevel(graph)
    
    TD = []
    for i, j in zip(A.nonzero()[0], A.nonzero()[1]):
        TD.append(s[i].item() - s[j].item())
    TDmatrix = coo_matrix((TD, (A.nonzero()[0], A.nonzero()[1])))
    return TDmatrix.toarray()


def HierarchicalLD(graph):
    A = nx.adjacency_matrix(graph).transpose()
    s = HierarchicalLevel(graph)
    
    TD = []
    for i, j in zip(A.nonzero()[0], A.nonzero()[1]):
        TD.append(s[i].item() - s[j].item())
    TDmatrix = coo_matrix((TD, (A.nonzero()[0], A.nonzero()[1])))
    return s, TDmatrix.toarray()


def HierarchicalCoherence(graph):
    A = nx.adjacency_matrix(graph).transpose()
    s = HierarchicalLevel(graph)
    TD = []
    for i, j in zip(A.nonzero()[0], A.nonzero()[1]):
        TD.append(s[i].item() - s[j].item())
    return HierarchicalDifferences(graph), mean(TD), pstdev(TD)



# returns a measure of equitable controllability over the full graph/network
def DemocracyCoefficient(graph):
    A = nx.adjacency_matrix(graph).transpose()
    s = HierarchicalLevel(graph)
    TD = []
    for i, j in zip(A.nonzero()[0], A.nonzero()[1]):
        TD.append(s[i].item() - s[j].item())
    return 1 - mean(TD)


#Do we just want the mean of the hierarchical differences over a node or do we also want the set of differences

# This is applying the democracy coefficient to a single node and is a measure of that nodes influence
# If the function returns masked that means the influence centrality of the node is 1, i.e basal node 
# The node variable call in the function is the index position/networkx numrical label of the nodes
# If node label/node index is out of bounds this also means the influence centrality is zero
def InfluenceCentrality(graph, node):
    masked = np.ma.masked_equal(HierarchicalDifferences(graph), 0)
    IC = masked[node].mean()
    return 1 - IC

     