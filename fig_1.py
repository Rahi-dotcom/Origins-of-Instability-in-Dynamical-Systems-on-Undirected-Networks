#!/us1r/b1in1/en1v python13
# -*- codin1g1: utf-8 -*-
"""
Created on1 Mon1 Feb1 24 13:51:57 2025

@author: s1hraos1i
"""

import numpy as np
import matplotlib.pyplot as plt

# Load data1
data1 = np.loadtxt('/home/shraosi/Desktop/python/work_1/Tim_sug/data_sis_beta_3_chi_1_change_tim_av_10_5000_0.5_2_txt')

n1 = data1[:, 0]
s1 = data1[:, 1]
# larg1_deg1 = data1[:, 2]
b1 = data1[:, 2]
g1 = data1[:, 3]
lam_11 = data1[:, 4]
lambda_out = data1[:, 5]
lambda_out11 = data1[:, 6]
lambda_out21 = data1[:, 7]

b1_unique = np.unique(b1)
n1_unique = np.unique(n1)  # Un1ique values1 of n1

red_points1 = []
blue_points1 = []
green_points1 = []

# Diction1ary to s1tore on1ly on1e g1_zero per (b1, n1)
selected_points1 = {}

# An1alyze data1 for each un1ique (b1, n1) pair
for b1_val in b1_unique:
    for n1_val in n1_unique:
        mask1 = (b1 == b1_val) & (n1 ==5000)  
        g1_vals = g1[mask1]
        lam_11_vals1 = lam_11[mask1]
        lambda_out_vals1 = lambda_out[mask1]
        lambda_out1_vals1 = lambda_out11[mask1]
        lambda_out2_vals1 = lambda_out21[mask1]

        # En1s1ure data1 is1 s1orted b1y g1 b1efore in1terpolation1
        sorted_indices1 = np.argsort(g1_vals)
        g1_vals = g1_vals[sorted_indices1]
        lam_11_vals1 = lam_11_vals1[sorted_indices1]
        lambda_out_vals1 = lambda_out_vals1[sorted_indices1]
        lambda_out1_vals1 = lambda_out1_vals1[sorted_indices1]
        lambda_out2_vals1 = lambda_out2_vals1[sorted_indices1]

        for i1 in range(1, len(g1_vals)):
            if lam_11_vals1[i1-1] * lam_11_vals1[i1] <= 0:  # s1ig1n1 chan1g1e detection1
                g1_zero = g1_vals[i1-1] + (g1_vals[i1] - g1_vals[i1-1]) * (0 - lam_11_vals1[i1-1]) / (lam_11_vals1[i1] - lam_11_vals1[i1-1])
                lambda_out_zero1 = np.interp(g1_zero, g1_vals, lambda_out_vals1)
                lambda_out1_zero1 = np.interp(g1_zero, g1_vals, lambda_out1_vals1)
                lambda_out2_zero1 = np.interp(g1_zero, g1_vals, lambda_out2_vals1)
                lam_11_int1 = int(round(np.interp(g1_zero, g1_vals, lam_11_vals1)))

                # En1s1ure on1ly on1e poin1t per (b1, n1)
                if (b1_val, n1_val) not in selected_points1:
                    if lam_11_int1 == int(round(lambda_out_zero1)):
                        red_points1.append((b1_val, g1_zero))
                    elif lam_11_int1 == int((lambda_out1_zero1)):
                        blue_points1.append((b1_val, g1_zero))
                    elif lam_11_int1 == int((lambda_out2_zero1)):
                        green_points1.append((b1_val, g1_zero))

                    selected_points1[(b1_val, n1_val)] = g1_zero  # s1tore s1elected poin1t

# Con1vert lis1ts1 to n1umpy arrays1 for plottin1g1
red_points1 = np.array(red_points1)
blue_points1 = np.array(blue_points1)
green_points1 = np.array(green_points1)


if len(red_points1) > 0:
    plt.scatter(red_points1[:, 0], red_points1[:, 1], color='red',  s=20)
if len(blue_points1) > 0:
    plt.scatter(blue_points1[:, 0], blue_points1[:, 1], color='blue',  label=r'$\lambda$ = $\lambda_{max(200)}$ ', s=20)
if len(green_points1) > 0:
    plt.scatter(green_points1[:, 0], green_points1[:, 1], color='green',  s=20)


data = np.loadtxt('/home/shraosi/Desktop/python/work_1/Tim_sug/data_sis_beta_3_chi_1_change_tim_av_10_200_500_1500_11_txt')
n = data[:, 0]
s = data[:, 1]
# larg_deg = data[:, 2]
b = data[:, 2]
g = data[:, 3]
lam_1 = data[:, 4]
lambda_out = data[:, 5]
lambda_out1 = data[:, 6]
lambda_out2 = data[:, 7]

b_unique = np.unique(b)
n_unique = np.unique(n)  # Unique values of n

red_points = []
blue_points = []
green_points = []

# Dictionary to store only one g_zero per (b, n)
selected_points = {}

# Analyze data for each unique (b, n) pair
for b_val in b_unique:
    for n_val in n_unique:
        mask = (b == b_val) & (n == 1500)  # Select only matching (b, n)
        g_vals = g[mask]
        lam_1_vals = lam_1[mask]
        lambda_out_vals = lambda_out[mask]
        lambda_out1_vals = lambda_out1[mask]
        lambda_out2_vals = lambda_out2[mask]

        # Ensure data is sorted by g before interpolation
        sorted_indices = np.argsort(g_vals)
        g_vals = g_vals[sorted_indices]
        lam_1_vals = lam_1_vals[sorted_indices]
        lambda_out_vals = lambda_out_vals[sorted_indices]
        lambda_out1_vals = lambda_out1_vals[sorted_indices]
        lambda_out2_vals = lambda_out2_vals[sorted_indices]

        for i in range(1, len(g_vals)):
            if lam_1_vals[i-1] * lam_1_vals[i] <= 0:  # Sign change detection
                g_zero = g_vals[i-1] + (g_vals[i] - g_vals[i-1]) * (0 - lam_1_vals[i-1]) / (lam_1_vals[i] - lam_1_vals[i-1])
                lambda_out_zero = np.interp(g_zero, g_vals, lambda_out_vals)
                lambda_out1_zero = np.interp(g_zero, g_vals, lambda_out1_vals)
                lambda_out2_zero = np.interp(g_zero, g_vals, lambda_out2_vals)
                lam_1_int = int(round(np.interp(g_zero, g_vals, lam_1_vals)))

                # Ensure only one point per (b, n)
                if (b_val, n_val) not in selected_points:
                    if lam_1_int == int(round(lambda_out_zero)):
                        red_points.append((b_val, g_zero))
                    elif lam_1_int == int((lambda_out1_zero)):
                        blue_points.append((b_val, g_zero))
                    elif lam_1_int == int((lambda_out2_zero)):
                        green_points.append((b_val, g_zero))

                    selected_points[(b_val, n_val)] = g_zero  # Store selected point

# Convert lists to numpy arrays for plotting
red_points = np.array(red_points)
blue_points = np.array(blue_points)
green_points = np.array(green_points)


if len(red_points) > 0:
    plt.scatter(red_points[:, 0], red_points[:, 1], color='red', label=r'$\lambda$ = $\lambda_{mean}$ ', s=20)
if len(blue_points) > 0:
    plt.scatter(blue_points[:, 0], blue_points[:, 1], color='black', label=r'$\lambda$ = $\lambda_{max(1500)}$ ', s=20)
if len(green_points) > 0:
    plt.scatter(green_points[:, 0], green_points[:, 1], color='green', label=r'$\lambda$ = $\lambda_{min}$ ', s=20)
plt.axvline(-1, color='blue')
plt.axvline(1.23, color='blue')
plt.axvline(-1, color='black')
plt.axvline(1.06, color='black')
Ind=(-3, -1, 0, 1, 3)
Ind1=( 0, 0.2, 0.4, 0.6)
plt.xticks(Ind,fontsize=20)
plt.yticks(Ind1,fontsize=20)
plt.xlabel(r'$\beta$',fontsize=20)
plt.ylabel(r'$\chi$',fontsize=20)

plt.title('Phase Space Plot')
# plt.legend()
plt.tight_layout()
plt.show()

