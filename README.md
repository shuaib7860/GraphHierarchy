# Graph Hierarchy

**Citation info**:  Moutsinas, G., Shuaib, C., Guo, W., & Jarvis, S. (2019), Graph hierarchy and spread of infections, arXiv preprint arXiv:1908.04358

## Overview
Trophic levels and hence trophic coherence can be defined only on networks with well defined sources and so trophic analysis of networks had been restricted to the ecological domain, until now. Graph Hierarchy is a python package that allows for analysis of all network structures via the trophic levels and coherence approach.  Trophic coherence, a measure of a network’s hierarchical organisation, has been shown to be linked to a network’s structural and dynamical aspects. In GraphHierarchy we have developed the python code which implements the mathematical generalisation of the trophic coherence theory to all networks. See citation paper for more details. 

## Features
  - Calculate hierarchical levels for any graph/network
  - Calculate the matrix of hierarchical differences for any graph/network
  - Calculate the mean and standard deviation of the distribution of hierarchical differences 
  - Calculate the democracy coefficient for any graph
  - Calculate the influence centrality for any node in a graph or for every single node in the graph in a single instance of code.


### Installation

The dependencies of the GraphHierarchy module are numpy, scipy, networkx and the statistics module. GraphHierarchy can be installed via the pip command.

```sh
pip install GraphHierarchy
```
And of course GraphHierarchy itself is open source with a [public repository] on GitHub. GraphHierarchy requires Python 3.5 or above to run. 

## How to use GraphHierarchy

Let's import the module after you have installed it, see installation instructions above for more details on how to install the GraphHierarchy module.

Firstly we create an instance of a networkx graph which will be the network we analyse with our GraphHierarchical functions.

```sh
import networkx as nx
import GraphHierarchy as gh
graph = nx.gnr_graph(20, 0.4)
nx.draw_networkx(graph)
```
The fourth line of code in the above script is a function call to visualise the graph. To calculate the hierarchical levels all one needs to write is:

```sh
gh.hierarchical_levels(graph, None)
```
This returns an n-dimensional array of the hierarchical levels in the ordering and with the labelling of graph.nodes(). If one wants to calculate the hierarchical differences and associated mean and standard deviation of this distribution:

```sh
gh.hierarchical_differences(graph, None)
```

This returns a three element tuple, the first element is an nxm dimensional array representing the hierarchical differences as an adjacency matrix for all the non-zero values. The second and third element of the tuple are the mean and the standard deviation of this distribution of differences respectively. 

If one requires the hierarchical levels, hierarchical differences and associated mean and standard deviation in one go:

```sh
gh.hierarchical_coherence(graph, None)
```

This returns a four element tuple, the first element is an n dimensional array representing the hierarchical levels, the second element is an nxm dimensional array representing the hierarchical differences as an adjacency matrix for all the non-zero values. The third and fourth element of the tuple are the mean and the standrard deviation of this distribution of differences respectively. 

We can also calculate the democracy coefficient for a graph:

```sh
gh.democracy_coefficient(graph, None)
```

This returns a single value, which is a topological metric. We can also work out the influence centrality for the node of a graph:

```sh
gh.influence_centrality(graph, None, node)
```
this returns a single numerical value. We can also work out the influence centrality for every single node in the graph:

```sh
gh.total_influence_centrality(graph, None)
```

The return of this function is a list of numerical values which are the influence centralities for nodes ordered according to graph.nodes(). This exhausts all the functions that are currently available in the module but we hope to add some more in newer versions of the module. 


# More Info
Fore more details on the mathematics behind these graph metrics see this [article] on arxiv. 

License
----

MIT

[article]: <https://arxiv.org/abs/1908.04358>
[public repository]: <https://github.com/shuaib7860/GraphHierarchy>
