#!/bin/python
import numpy as np
import matplotlib.pylab as plt
import seaborn as sbn
import sys
import os
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import RMSD
sbn.set(font_scale=2,style="white")

# Produces conformation plots (syn vs anti) for MIF 180
# for all stages in a specified run (given as an argument). Run from
# base directory
# Note this will only work if minimal coordinate saving was set to
# True in sim.cfg

run = 1
#run = int(sys.argv[1])

# Define atom selection
selection = 'resname LIG and (not name H* OAA CAS CAD CAF CAG CAP CAE)'

# define function to extract data and create plots

def plot(subfig ,leg, stage):

    #Define reference
    ref = mda.Universe(f'{leg}/run00{run}/{stage}/input/SYSTEM.pdb')

    #Define lambda values
    lam_folders = [x for x in os.listdir(f'{leg}/run00{run}/{stage}/output') if 'lambda' in x]
    lam_folders.sort()

    #Plot rmsd for all windows
    no_windows = len(lam_folders)

    # Add subtitle
    subfig.suptitle(f'{leg} {stage}')

    for i in range(no_windows):
        
        #calculate RMSD for window
        mobile = mda.Universe(f'{leg}/run00{run}/{stage}/input/SYSTEM.pdb', 
                             f'{leg}/run00{run}/{stage}/output/{lam_folders[i]}/traj000000001.dcd')
        R = RMSD(mobile, ref, select=selection, groupselections=[selection])
        R.run()
        rmsd = R.results.rmsd.T
        time = rmsd[1]
        
        #plot RMSD for window
        ax = subfig.add_subplot(no_windows, 1, i+1)
        if True in (rmsd[3] > 0.7):
            ax.plot(time, rmsd[3], 'r-',label=lam_folders[i])
            print(f'{lam_folders[i]} ({leg} {stage}) window samples opposite conformation')
        else:
            ax.plot(time, rmsd[3], 'k-',label=lam_folders[i])
        
        ax.legend(loc="best")
        
        if i!= no_windows - 1:
            x_axis = ax.get_xaxis()
            x_axis.set_visible(False)

        if lam_folders[i] == lam_folders[-1]:
            ax.set_xlabel("time (ps)")
        
        if leg == "bound" and stage == "discharge":
             ax.set_ylabel(r"RMSD ($\AA$)")

# plot conformation plots

fig, ax1 = plt.subplots(1,1,figsize=(32,78))
# ax1.text(0.5, 0.05, 'Time (ps)', ha='center',fontsize = 'x-large')
# ax1.text(0.05, 0.5, 'RMSD ($\AA$)', va='center', rotation='vertical', fontsize='x-large')

subfigs = fig.subfigures(1,4)

i = 0
for leg in ["bound", "free"]:
    for stage in ["discharge", "vanish"]:
        plot(subfigs[i], leg=leg, stage=stage)
        i += 1

fig.savefig(f'conformation-plots-run-{run}.png', pad_inches=3)