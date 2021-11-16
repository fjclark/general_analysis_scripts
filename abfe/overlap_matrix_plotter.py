#!/bin/python
import numpy as np
import matplotlib.pylab as plt
import seaborn as sbn
import sys
sbn.set(font_scale=2,style="white")

#read the data file
f = open(sys.argv[1], 'r')
lines = f.readlines()

# Calculate no lamda values
line_idx_overlap_matrix = 0
line_idx_dg = 0

for i,l in enumerate(lines):
    if l.startswith('#Overlap matrix'):
        line_idx_overlap_matrix = i
    if l.startswith('#DG from neighbouring lambda in kcal/mol'):
        line_idx_dg = i

no_lam_vals = line_idx_dg - line_idx_overlap_matrix -1

#extract data
matrix = []
pmf_ti = []
pmf_mbar = []
for i,l in enumerate(lines):
    if l.startswith('#Overlap matrix'):
        for int in range(i+1, i+no_lam_vals+1):
            matrix.append(np.array(lines[int].strip().split(' ')).astype('float'))
    elif l.startswith('#PMF from MBAR in kcal/mol'):
        for int in range(i+1, i+no_lam_vals+1):
            pmf_mbar.append(np.array(lines[int].strip().split(' ')).astype('float'))
    elif l.startswith('#PMF from TI in kcal/mol'):
        for int in range(i+1, i+no_lam_vals+1):
            pmf_ti.append(np.array(lines[int].strip().split(' ')).astype('float'))

matrix = np.array(matrix)
pmf_ti = np.array(pmf_ti)
pmf_mbar = np.array(pmf_mbar)

#plot pmf
plt.figure(figsize = (10,7))
plt.plot(pmf_ti[:,0],pmf_ti[:,1], lw=0.5,marker='o', label = 'TI')
plt.errorbar(pmf_mbar[:,0],pmf_mbar[:,1],yerr=pmf_mbar[:,2], elinewidth=2, lw=0.5,marker='o',label = 'MBAR')
plt.xlabel(r'$\lambda$')
plt.ylabel(r'$\Delta$ G in [kcal/mol]')
sbn.despine()
plt.legend()
plt.savefig('pmf_plot.png',dpi=400, bbox_inches = "tight")

#plot overlap matrix
plt.figure(figsize = (10,7))
sbn.heatmap(matrix, square=True, linewidths=.5).figure.savefig('overlap_plot.png',dpi=400, bbox_inches = "tight")
