# Synchrony Data Project - Lybrand Lab
## Logic
### Clustering Coefficient
Let $c_{ij}$ represent an existing connection (i.e. Synchronization Index above average) between node i and node j.
The following set, $N_i$, then, expressed all nodes, $n$, that are connected to node i, where C represents a set of all existing connections within the network:\\
$N_i = \{n_j:c_{ij} \in C \}$

For convenience, $|N_i| = l_i$. If node i is connected to 3 nodes, $l_i = 3$.\\
For node i, a connection $c_{ij}$ is referred to as a direct connection, while connection $c_{jk}$, where both $n_j$ and $n_k$ are elements of $N_i$, is an indirect connection.

The clustering coefficient of node i is calculated by observing the indirect connections, defined by the following expression:\\
$CC_i = \frac {|\{c_{jk} : n_j \in N_i, n_k \in N_i\}} {l_i (l_i-1)}$

Thus, a high Clustering Coefficient suggests that the nodes connected to node i are all strongly connected to each other.
For a node with a clustering coefficient of 1, every possible indirect exists. We call a node "significant" if its Clustering Coefficient is above average for its phase. A node is "perfect" if its Clustering coefficient is equal to 1. By the range of a Clustering Coefficient, all perfect nodes are significant, but not all significant nodes are necessarily perfect.

### Network
A perfect node's network is defined as the complete set of channels that connect to the perfect node. For example, if channel 12 is a perfect node and is connected to channels 13, 22, and 42, the set of nodes {12, 13, 22, 42} would be considered channel 12's network. By definition, all elements of a perfect node's network must be connected to all other elements of that network. Now, suppose channel 13 is also a perfect network. All of channel 13's connections must then be connected and all of channel 12's connections (that is all elements of channel 12's network) must be connected to channel 13. Thus, by definition of the Clustering Coefficient, we conclude that, if 2 or more perfect nodes are connected, their network must be identical. We record each perfect node's network as a distinct network. In the image below, we notice 2 distinct networks, one (red) centered around channels 12 and 13, and one centered around channel 32:

For the sake of consistency with contemporary Graph Theory conventions, we refer to these networks as Small-World Networks (SWNs).
