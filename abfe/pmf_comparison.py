#!/bin/python

def plot(datafile, ax ,leg, stage):

    #read the data file
    with open(datafile,'r') as f:
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
    ax.errorbar(pmf_mbar[:,0],pmf_mbar[:,1],yerr=pmf_mbar[:,2], elinewidth=2, lw=0.5,marker='o',label = f'run {run}')
    ax.set_xlabel(r'$\lambda$')
    ax.set_ylabel(r'$\Delta$ G in [kcal/mol]')
    sbn.despine()

    
    ax.set_title(f'{leg} {stage}')

if __name__ == "__main__":

    import numpy as np
    import matplotlib.pylab as plt
    import seaborn as sbn
    import sys
    sbn.set(font_scale=2,style="white")

    # Produces overlap matrices figure for all stages in a 
    # specified run (given as an argument)

    runs = int(sys.argv[1])
    #runs = 3

    #plot PMFs 

    fig, axs = plt.subplots(2,2, figsize=(12,12))

    for run in range(1, runs+1):
        #plot bound on top row
        for i, stage in enumerate(['discharge','vanish']):
            ax = axs[0,i]
            plot(f'bound/run00{str(run)}/{stage}/output/freenrg-MBAR-p-83-overlap.dat', ax, 'bound', stage)
        # plot free on bottom row
        for i, stage in enumerate(['discharge','vanish']):
            ax = axs[1,i]
            plot(f'free/run00{str(run)}/{stage}/output/freenrg-MBAR-p-83-overlap.dat', ax, 'free', stage)
            #set legend for last PMF plot
            if run == runs and stage == 'vanish':
                ax.legend()
            
        fig.tight_layout()
        fig.savefig('pmf_comparison.png')