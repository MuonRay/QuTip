# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 01:03:00 2022

@author: cosmi
"""
import networkx as nx
import numpy as np
import scipy
import scipy.spatial
import matplotlib.pyplot as plt

from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
n = 4
eta = 4

G = nx.scale_free_graph(n) #obtain a directed graph
#remove self loops
G.remove_edges_from(nx.selfloop_edges(G))
#remove parallel links

A = nx.adj_matrix(G).todense().T
k_in = np.zeros(G.number_of_nodes())
            
#weighted laplacian matrix
L = np.diag(k_in) - A
print(L)



