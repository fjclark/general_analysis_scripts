
#!/bin/python
import numpy as np
import matplotlib.pylab as plt
import seaborn as sbn
import sys
import os
import re
sbn.set(font_scale=2,style="white")

# Produces overlap matrices figure for all stages in a 
# specified run (given as an argument)

# 01 to allow for up to 99 runs

run_no = '0'

# input desired run number as a string
input_no = sys.argv[1]

if int(input_no) < 10:
    run_no += input_no
else:
    run_no = input_no

# Get folder name
dir_names = os.listdir("bound")
r = re.compile(f"run0{run_no}")
dir_name = list(filter(r.match, dir_names))[0]

# define function to extract data and plot PMFs/ overlap matrices

def plot(datafile, ax ,leg, stage, pmf_or_overlap="overlap"):

    #read the data file
    f = open(datafile,'r')
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

    if pmf_or_overlap == 'pmf':
        #plot pmf
        ax.plot(pmf_ti[:,0],pmf_ti[:,1], lw=0.5,marker='o', label = 'TI')
        ax.errorbar(pmf_mbar[:,0],pmf_mbar[:,1],yerr=pmf_mbar[:,2], elinewidth=2, lw=0.5,marker='o',label = 'MBAR')
        ax.set_xlabel(r'$\lambda$')
        ax.set_ylabel(r'$\Delta$ G in [kcal/mol]')
        sbn.despine()

    elif pmf_or_overlap == 'overlap':
        #plot overlap matrix
        sbn.heatmap(matrix, ax=ax, square=True, linewidths=.5).figure
    
    ax.set_title(f'{leg} {stage}')

#plot PMFs and Overlap

for plot_type in ['pmf', 'overlap']:
    fig = plt.figure(figsize=(12,12))
    #plot bound on top row
    for i, stage in enumerate(['discharge','vanish']):
        ax = fig.add_subplot(2,2,i+1)
        plot(f'bound/{dir_name}/{stage}/output/freenrg-MBAR-p-83-overlap.dat', ax, 'bound', stage, pmf_or_overlap=plot_type)
    # plot free on bottom row
    for i, stage in enumerate(['discharge','vanish']):
        ax = fig.add_subplot(2,2,i+3)
        plot(f'free/{dir_name}/{stage}/output/freenrg-MBAR-p-83-overlap.dat', ax, 'free', stage, pmf_or_overlap=plot_type)
        #set legend for last PMF plot
        if plot_type == 'pmf' and stage == 'vanish':
            ax.legend()
        
    fig.tight_layout()
    fig.savefig(f'convergence-run-{run_no}_individual-{plot_type}.png')
