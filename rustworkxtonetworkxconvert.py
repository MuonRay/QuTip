# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 21:47:19 2022

@author: cosmi
"""

import networkx as nx
import rustworkx as rx
from pylab import show

from rustworkx.visualization import mpl_draw

global g, nextg

g = rx.PyDiGraph(check_cycle=True)

# Each time add node is called, it returns a new node index
a = g.add_node("A")
b = g.add_node("B")
c = g.add_node("C")

# add_edges_from takes tuples of node indices and weights,
# and returns edge indices
g.add_edges_from([(a, b, 1.5), (a, c, 5.0), (b, c, 2.5)])

# Returns the path A -> B -> C
rx.dijkstra_shortest_paths(g, a, c, weight_fn=float)



def convert_rustworkx_to_networkx(g):
    """Convert a rustworkx PyGraph or PyDiGraph to a networkx graph."""
    edge_list = [(
        g[x[0]], g[x[1]],
        {'weight': x[2]}) for x in g.weighted_edge_list()]

    if isinstance(g, rx.PyGraph):
        if g.multigraph:
            return nx.MultiGraph(edge_list)
        else:
            return nx.g(edge_list)
        nx.draw(g)
        show()
    else:
        if g.multigraph:
            return nx.MultiDiGraph(edge_list)
        else:
            return nx.DiGraph(edge_list)


convert_rustworkx_to_networkx(g)
print(g)

    
