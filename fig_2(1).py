#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  5 14:22:42 2025

@author: shraosi
"""

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from scipy.sparse.linalg import eigsh

# Parameters
N = 200  # Number of nodes
m = 5  # Number of edges to attach per new node in BA model
mu, nu, rho = 1, -1, 0
chi_values = [-0.1, 0.01, 0.3]


# G = nx.barabasi_albert_graph(N, m)
# degrees = np.array([d for _, d in G.degree()])
B = np.loadtxt('/home/shraosi/Desktop/python/work_1/Tim_sug/BA_m_5_n_200.txt', float)
G = nx.Graph(B)
A = nx.to_numpy_array(G)
degrees = np.sum(A, axis=0)
theo = degrees**(nu+1)
theo /=  np.max(theo)

beta = 0  
J = np.zeros((N, N))


# for i in range(N):
#     J[i, i] = -beta - chi_values[1] * degrees[i]**mu  # Diagonal elements
#     for j in range(N):
#         # if A[i, j] == 1:
#             J[i, j] = degrees[i]**nu * (degrees[j]**rho)*A[i,j]  # Off-diagonal elements


# eigvals, eigvecs = eigsh(J, k=1, which='LA')  # Largest eigenvalue
# v = np.abs(eigvecs[:, 0])  # Take absolute values of eigenvector components
# v_th=[]

fig, axs = plt.subplots(1, 3, figsize=(12, 4), sharey=True)

for i, chi in enumerate(chi_values):
    J = np.zeros((N, N))
    for j in range(N):
        
        for k in range(N):
            
            J[j,k] = A[j,k]*degrees[j]**(nu)*degrees[k]**(rho)
            J[j, j] = -beta - chi * degrees[j]**mu
    
    eigvals, eigvecs = eigsh(J, k=1, which='LA')
    v = np.abs(eigvecs[:, 0])
    


    axs[i].scatter(degrees, v, c=v, cmap='coolwarm', alpha=0.7)
    axs[1].plot(np.sort(degrees), 0.07*np.sort(theo), color='red')#, linestyle='dashed', label=r'$v_i \sim k_i^{\nu+1}$')
    Ind=( 0, 20,40,60)
    Ind1=( 0, 0.2, 0.4, 0.6,0.8,1)
    axs[i].set_xticks(Ind)
    axs[i].set_xticklabels(Ind, fontsize=20)

    axs[i].set_yticks(Ind1)
    axs[i].set_yticklabels(Ind1, fontsize=20)
    # axs[i].set_title(f'χ = {chi}',fontsize=20)
    axs[i].set_xlabel(r'$d$',fontsize=20)
# v_th.append(v_i)

axs[0].set_ylabel(r'$v$',fontsize=20)
plt.tight_layout()
plt.show()
