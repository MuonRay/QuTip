# QuTip Quantum Many-Body and Network Simulations
a collection of codes for use in simulating quantum systems, with particular attention of generating Wigner function, Density matrix and other plots to monitor the number of states as they grow in complexity and evolve over time. 

Of particular interest is to model many-body networks as a cavity with a designated input and output class node(s) to effectively do benchmarking of quantum networks using trace distance, fidelity measures, etc. 

Also a significant portion of code used here is tailored to work with networkx and rustworkx for directional acyclic graphs that can be used to model quantum networks as networks of Kuramoto-like harmonic oscillators wherein the synchronization IS entanglement.
