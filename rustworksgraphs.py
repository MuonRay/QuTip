# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 21:16:23 2022

@author: cosmi

Rustworkx uses some of the similarities of networkx but gives the necessary criterion of acylic 
definition of graphs as quantum netowkrs differ from their classical counterparts
as having a clear "arrow of time" in terms of phenominon such as measurements, updates,
collapse of superposition states, wavefunction time evolution and so forth.
Rustworkx evolved from retworkx used by QISKIT for their operations of 
gate baed quantum computers using graph theory to internalised the gate operations
"""

import rustworkx as rx
from rustworkx.visualization import mpl_draw
#acyclic path (quantum)
path_graph = rx.generators.directed_path_graph(5)
mpl_draw(path_graph)


#cyclic path (classical)
#cycle_graph = rx.generators.directed_cycle_graph(5)
#mpl_draw(cycle_graph)

dag = rx.PyDiGraph()
# Input nodes:
in_nodes = dag.add_nodes_from(["q_0", "q_1", "c_0", "c_1"])
# Output nodes
out_nodes = dag.add_nodes_from(["q_0", "q_1", "c_0", "c_1"])
# Equivalency matrix
translation_matrix = {"h": ["rz(pi/2)", "sx", "rz(pi/2)"]}
# Insructions natively supported on target QPU
hardware_instructions = {"measure", "cx", "sx", "rz", "x"}

# Iterate over instructions in order and replace gates outside of native
# instruction set with a subcircuit from the translation matrix
for gate_index in rx.topological_sort(dag):
    if gate_index not in in_nodes and gate_index not in out_nodes:
        if dag[gate_index] not in hardware_instructions:
            edge_val = dag.out_edges(gate_index)[0][2]
            equivalent_subcircuit = rx.PyDiGraph()
            count = 0
            for node in translation_matrix[dag[gate_index]]:
                if count == 0:
                    equivalent_subcircuit.add_node(node)
                else:
                    equivalent_subcircuit.add_child(count - 1, node, edge_val)
                count += 1

            def map_fn(source, target, weight):
                if source == gate_index:
                    return len(equivalent_subcircuit) - 1
                else:
                    return 0

            dag.substitute_node_with_subgraph(
                gate_index,
                equivalent_subcircuit,
                map_fn
            )



bit_nodes = {"q_0", "q_1", "c_0", "c_1"}

def filter_fn(node):
    # Don't collect input or output nodes
    if node in bit_nodes:
        return False
    # Don't include 2 qubit gates
    if node == "cx":
        return False
    # Ignore non-unitary operations
    if node == "measure":
        return False
    return True

runs = rx.collect_runs(dag, filter_fn)
print(runs)