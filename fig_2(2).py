#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 17 11:13:32 2025

@author: shraosi
"""
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from scipy.sparse.linalg import eigsh

# Parameters
N = 200  # Number of nodes
m = 5    # Number of edges to attach per new node in BA model
nu = -1
rho = 0
mu = 1
beta = 0
chi_values = [-0.1, 0.01, 0.3]  # Different values of chi

# Load adjacency matrix
B = np.loadtxt('/home/shraosi/Desktop/python/work_1/Tim_sug/BA_m_5_n_200.txt', float)
G = nx.Graph(B)
A = nx.to_numpy_array(G)
degrees = np.sum(A, axis=0)

# Compute 50th percentile threshold (median degree)
degree_threshold = np.percentile(degrees, 50)

# Filter nodes with degrees above threshold
high_degree_nodes = [node for node in G.nodes() if degrees[node] > degree_threshold]
H = G.subgraph(high_degree_nodes)  # Induced subgraph

# Compute new positions for the filtered network
pos = nx.spring_layout(H)

# Plot for different chi values
fig, axes = plt.subplots(1, 3, figsize=(12, 4))

for i, chi in enumerate(chi_values):
    J = np.zeros((N, N))
    for j in range(N):
        for k in range(N):       
            J[j, k] = (degrees[j]**nu) * (degrees[k]**rho) * A[j, k] 
        J[j, j] = -beta - chi * degrees[j]**mu

    # Compute largest eigenvector
    eigvals, eigvecs = eigsh(J, k=1, which='LA')  
    largest_eigenvector = eigvecs[:, 0]

    # Normalize colors and node sizes
    node_colors = {node: largest_eigenvector[node] for node in high_degree_nodes}
    node_sizes = {node: (degrees[node] + 1) * 15 for node in high_degree_nodes}

    # Sort nodes by eigenvector values for layering effect
    sorted_nodes = sorted(high_degree_nodes, key=lambda n: node_colors[n])

    ax = axes[i]
    
    # Draw low-eigenvector nodes first
    nx.draw_networkx_nodes(
        H, pos, 
        nodelist=sorted_nodes, 
        node_color=[node_colors[n] for n in sorted_nodes], 
        cmap='coolwarm', 
        node_size=[node_sizes[n] for n in sorted_nodes], 
        edgecolors='black', alpha=0.7, ax=ax
    )
    
    # Draw edges
    nx.draw_networkx_edges(H, pos, edge_color='black', alpha=0.3, ax=ax)

    # ax.set_title(f'chi = {chi}')

plt.tight_layout()
plt.show()
