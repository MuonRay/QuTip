# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 02:01:45 2022

@author: cosmi
Takes the adjacency matrix of a UCG and the indices of the
starting nodes, returns the dictionary of a DAG."""

def ucg2dag(adj_matrix, starts):
    adj_list = [
        [target for target, is_connected in enumerate(row) if is_connected]
            for row in adj_matrix
    ]

    frontier = starts

    dag = [[] for _ in range(len(adj_list))]

    while frontier:
        for source in frontier:
            dag[source].extend(target for target in adj_list[source] if not target in starts)
        frontier = set(target 
            for source in frontier for target in adj_list[source] if not target in starts
        )
        starts.update(frontier)

    return dag