#!/usr/bin/env python
# coding: utf-8

import matplotlib.pyplot as plt
import os
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import RMSD

#Define reference and atom selection
ref = mda.Universe('../input/SYSTEM.pdb')
selection = 'resname 3TX and (not name H OAA CAS CAD CAF CAG CAP CAE)'

#Define lambda values
lam_folders = [x for x in os.listdir() if 'lambda' in x]
lam_folders.sort()

#Plot rmsd for all windows
fig = plt.figure(figsize=(6,24))
no_windows = len(lam_folders)

for i in range(no_windows):
    
    #calculate RMSD for window
    mobile = mda.Universe('../input/SYSTEM.pdb', lam_folders[i]+'/traj000000001.dcd')
    R = RMSD(mobile, ref, select=selection, groupselections=[selection])
    R.run()
    rmsd = R.results.rmsd.T
    time = rmsd[1]
    
    #plot RMSD for window
    ax = fig.add_subplot(no_windows, 1, i+1)
    if True in (rmsd[3] > 0.7):
        ax.plot(time, rmsd[3], 'r-',label=lam_folders[i])
        print(f'{lam_folders[i]} window samples syn conformation')
    else:
        ax.plot(time, rmsd[3], 'k-',label=lam_folders[i])
    ax.set_ylabel(r"RMSD ($\AA$)")
    ax.legend(loc="best")
    ax.set_xlabel("time (ps)")
    
    if i!= no_windows - 1:
        x_axis = ax.get_xaxis()
        x_axis.set_visible(False)
    
fig.savefig("conformation_rmsds.png")
