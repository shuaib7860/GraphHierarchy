import numpy as np 
import networkx as nx
from statistics import mean, pstdev
from scipy.sparse import coo_matrix, csr_matrix, diags
from scipy.sparse.linalg import lsqr