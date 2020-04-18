# Graph Hierarchy

**Citation info**:  Moutsinas, G., Shuaib, C., Guo, W., & Jarvis, S. (2019), Graph hierarchy and spread of infections, arXiv preprint arXiv:1908.04358

## Overview
Trophic levels and hence trophic coherence can be defined only on networks with well defined sources. Trophic coherence, a measure of a network’s hierarchical organisation, has been shown to be linked to a network’s structural and dynamical aspects. Thus trophic analysis of networks had been restricted to the ecological domain, until now. In GraphHierarchy we have created the python code which implements the mathematical generalisation of the trophic coherence theory to all networks. See citation info for more details. 

## Features
  - Calculate hierarchical levels for any graph
  - Calculate the weighted adjacency matrix of hierarchical differences for any graph
  - Calculate the mean and standard deviation of the distribution of hierarchical differences 
  - Calculate the democracy coefficient for any graph
  - Calculate the influence centrality for any node in a graph or for every single node in the graph.


### Installation

The dependencies of the GraphHierarchy module are numpy, scipy and networkx. GraphHierarchy can be installed via the pip command.

```sh
pip install GraphHierarchy
```
And of course GraphHierarchy itself is open source with a [public repository] on GitHub. GraphHierarchy requires Python 3.5 or above to run. 

## How to use GraphHierarchy

Let's import the module once you have installed it, see installation instructions above for more details on how to install the GraphHierarchy module.

Firstly we create an instance of a networkx graph which will be the network we analyse utilising our GraphHierarchical functions.

```sh
import networkx as nx
import GraphHierarchy as gh
graph = nx.gnr_graph(20, 0.4)
nx.draw_networkx(graph)
```
The fourth line of code in the above script is a function call to visualise the graph. To calculate the hierarchical levels, there are two parameters required, a network object, and a weight parameter. In the case of unweighted network all one needs to write is:

```sh
gh.hierarchical_levels(graph, None)
```
This returns an $N \times 1$ dimensional array of the hierarchical levels indexed by the nodes in the ordering of graph.nodes. When there is a weight parameter, the corresponding attribute should be input as a string as the second parameter in the function call. If one wants to calculate the hierarchical differences:

```sh
gh.hierarchical_differences(graph, None)
```

This returns a weighted adjacency matrix, an $N \times N $ dimensional array with the weights representing the hierarchical differences. 

The function hierarchical coherence returns the hierarchical differences and associated mean and standard deviation of this distribution:

```sh
gh.hierarchical_coherence(graph, None)
```

This returns a three element tuple, the first element is the hierarchical differences adjacency matrix. The second and third elements of the tuple are the mean of this distribution of hierarchical differences and the standard deviation of this distribution of differences respectively. The standard deviation of the distribution is known as the hierarchical coherence and is an important metric which gives a measure of a network's organisation and structure. 

We can also calculate the democracy coefficient, a topological metric, for a graph:

```sh
gh.democracy_coefficient(graph, None)
```

This returns a single value. We can also work out the influence centrality for the node of a graph:

```sh
gh.node_influence_centrality(graph, None, node)
```
This returns a single numerical value. We can also work out the influence centrality for every single node in the graph:

```sh
gh.influence_centrality(graph, None)
```

The return of this function is an $N \times 1$ dimensional array indexed by the nodes with values of their respective influence centralities in the same order as graph.nodes(). A final function returns all the possible node, edge and graph metrics so far developed in this theory.

```sh
gh.hierarchical_metrics(graph, None)
```

This exhausts all the functions that are currently available in the module but we hope to add some more in newer versions of the module. 


# More Info
Fore more details on the mathematics behind these graph metrics see this [article] on arxiv. 

License
----

MIT

[article]: <https://arxiv.org/abs/1908.04358>
[public repository]: <https://github.com/shuaib7860/GraphHierarchy>
