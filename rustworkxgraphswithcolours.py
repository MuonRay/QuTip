# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 19:42:29 2022

@author: cosmi
"""

import rustworkx
from rustworkx.visualization import mpl_draw

graph = rustworkx.generators.grid_graph(5, 5)
btw = rustworkx.betweenness_centrality(graph)
# Color nodes in graph visualization with betweenness centrality
colors = []
for i in range(len(graph)):
    colors.append(btw[i])
mpl_draw(graph, node_color=colors)