#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 12 04:05:48 2025

@author: shraosi
"""

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import eig
from scipy.integrate import solve_ivp
import math

mu_net = 1
chi = 1
alpha = 1

n = 500
p = 0.02
G = nx.erdos_renyi_graph(n, p, directed=True)

A = nx.to_numpy_array(G, dtype=int)
G = nx.Graph(A)

# Extract the giant component
connected_components = nx.connected_components(G)
giant_component = max(connected_components, key=len)
G1 = G.subgraph(giant_component)

def calculate_mean_non_zero_off_diagonal(matrix):
    n = matrix.shape[0]
    non_zero_off_diagonal_terms = [matrix[i, j] for i in range(n)
                                   for j in range(n) if i != j and matrix[i, j] != 0]
    return np.mean(non_zero_off_diagonal_terms)

A = nx.to_numpy_array(G1)
n = len(A)
num_edges = G1.number_of_edges()
num_nodes = G1.number_of_nodes()
degree = np.sum(A, axis=0)
avg = np.mean(degree)

total_possible_edges = (num_nodes * (num_nodes - 1)) / 2
s = num_edges / total_possible_edges

def calculate_degree_sparsity(G):
    degrees = [deg for _, deg in G.degree()]
    mean_k = np.mean(degrees)
    mean_k_squared = np.mean([k**2 for k in degrees])
    sparsity = mean_k_squared / mean_k**2
    return sparsity

sparsity = calculate_degree_sparsity(G1)

data = np.arange(0.01, 0.5, 0.02)
lam_1 = []
lam_t = []
lam_P1 = []
lam_P2 = []
lam_11 = []
lam_t1 = []
lam_P1 = []
lam_P2 = []

# Dynamical equation
def dxdt(t, x, beta, A):
    return -x + beta * (1 - x) * (A @ x)

for m in data:
    beta = m
    print("Processing β =", beta)

    # === Analytical matrices (keep same) ===
    x_avg = 1 - alpha / (beta * avg)
    M = np.zeros([n, n], float)
    for i in range(n):
        for j in range(n):
            M[i, j] = beta * A[i, j]
        M[i, i] = -alpha

    values, vect = eig(M)
    x1 = values.real
    x2 = np.sort(x1)
    result = calculate_mean_non_zero_off_diagonal(M)
    lambda_out = -alpha + result * (n - 1) * s * sparsity
    lam_t.append(lambda_out)
    lam_1.append(x2[-1])

    M1 = np.zeros([n, n], float)
    for i in range(n):
        for j in range(n):
            M1[i, j] = (beta * alpha / (alpha + beta * x_avg * degree[i])) * A[i, j]
        M1[i, i] = -alpha - beta * x_avg * degree[i]

    values1, vect = eig(M1)
    x11 = values1.real
    x21 = np.sort(x11)
    result1 = calculate_mean_non_zero_off_diagonal(M1)
    lambda_out1 = -alpha - x_avg * beta * avg + result1 * (n - 1) * s * sparsity
    lambda_out11 = -alpha - beta * x_avg * max(degree)
    lambda_out21 = -alpha - beta * x_avg * min(degree)
    lam_t1.append(lambda_out1)
    lam_11.append(x21[-1])
    lam_P1.append(lambda_out11)
    lam_P2.append(lambda_out21)

    # === New Part: Numerical steady state and Jacobian ===
    # Integrate dynamics to steady state
    x0 = np.random.rand(n)
    sol = solve_ivp(dxdt, [0, 500], x0, args=(beta, A), rtol=1e-6, atol=1e-8)
    x_star = sol.y[:, -1]

    # Jacobian at steady state
    J = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                J[i, j] = -1 - beta * (A @ x_star)[i] + beta * (1 - x_star[i]) * A[i, j]
            else:
                J[i, j] = beta * (1 - x_star[i]) * A[i, j]

    eigvals, _ = eig(J)
    lam_numerical = np.sort(eigvals.real)[-1]
    lam_11[-1] = lam_numerical  # replace previous analytical eigenvalue

    # For lam_1 also use numerical steady-state Jacobian
    lam_1[-1] = lam_numerical

# === Plotting ===
plt.plot(data, lam_t1, 'r')
plt.plot(data, lam_t, 'r')
plt.plot(data, lam_P2, 'b')
plt.plot(data, lam_11, 'k.')
plt.plot(data, lam_1, 'k.')
plt.axhline(y=0, color='black')
plt.xlabel(r'$\beta$', fontsize=15)
plt.ylabel(r'$\lambda_{\max}$', fontsize=15)
plt.show()
