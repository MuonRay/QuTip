# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 19:55:24 2022

@author: cosmi
"""
import rustworkx
import numpy as np
import matplotlib.pyplot as plt

import networkx as nx

from rustworkx.visualization import mpl_draw


matrix = np.array([[0, -1, -1], [1, 0, -1], [1, 1, 0]], dtype=np.float64)
graph = rustworkx.PyDiGraph.from_adjacency_matrix(matrix)
mpl_draw(graph, with_labels=True, edge_labels=str)



# Plot the Adjacency and the Laplacian matrices corresponding to the ring graph
fig,axs = plt.subplots(1,3, figsize=(15,3))
im = axs[0].imshow(matrix)
plt.colorbar(im,ax = axs[0])
axs[0].set_title('Adjacency matrix')
im = axs[1].imshow(matrix[:20,:20])
plt.colorbar(im,ax = axs[1])
axs[1].set_title('Adjacency matrix (zoom)')
