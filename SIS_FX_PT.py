# -*- coding: utf-8 -*-
"""
Created on Thu Apr 10 00:28:27 2025

@author: Shraosi Dawn
"""

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import scipy.io

# N = 500  
# p = 0.02  
beta_values = np.concatenate((np.linspace(0, 0.3, 50), np.linspace(0.31,0.5,10)))    
T_max = 10000  
t_eval = np.linspace(0, T_max, 10000)

theo_avg = []
theo_star = []
# G = nx.erdos_renyi_graph(N, p)
# A = nx.to_numpy_array(G)
network_data = scipy.io.loadmat('C:/Users/Shraosi Dawn/Desktop/NET PHYSICS/python/work_1/Tim_sug/fraction_node/A_BA_avgdeg12.mat')
mat_file = network_data['A']
G = nx.Graph(mat_file)
A = nx.to_numpy_array(G)
A = np.array(A)

N = len(A)

degree=np.sum(A, axis=0)
avg=np.mean(degree)


def dxdt(t, x, beta, A):
    return -x + beta * (1 - x) * (A @ x)


x_star_values = []

for beta in beta_values:
    print(beta)
    
    x0 = np.random.rand(N)

    
    sol = solve_ivp(dxdt, [0, T_max], x0, args=(beta, A), t_eval=t_eval, rtol=1e-6, atol=1e-8)
    
    
    x_final = sol.y[:, -1]
    x_star = np.mean(x_final)
    x_star_values.append(x_star)
    
    x_avg= 1 - 1/(beta*avg)
    
    x_star = 1 -1/( 1+(beta-1/avg)*avg)
    theo_avg.append(x_avg)
    theo_star.append(x_star)

ax = plt.gca()  
for spine in ax.spines.values():
    spine.set_edgecolor('black')
    spine.set_linewidth(3)


ax.tick_params(direction='in', colors='black', width=1.5)
ind1=(0,0.1,0.2,0.3,0.4,0.5)
ind=(0,0.2,0.4,0.6,0.8,1)
plt.xticks(ind1,fontsize=20,fontweight='bold')
plt.yticks(ind,fontsize=20,fontweight='bold')


#plt.scatter(beta_values ,theo_avg, color = 'red')
plt.scatter(beta_values, theo_star,color= 'red',s=70)
plt.plot(beta_values, x_star_values,lw=5)
plt.xlabel(r'$\beta$', fontsize=25,fontweight='bold')
plt.ylabel(r' $\langle x^*\rangle$', fontsize=25,fontweight='bold')
# plt.title(r'$\beta$ vs Mean Fixed Point $x^*$', fontsize=16)
#plt.grid(True)
plt.tight_layout()
plt.show()
