from cProfile import label


def avg_water(run_folder, lam_val, percent_frames, index, length):
    '''Returns the average number of water molecules within a given length of
    a given atom index for a given system (pdb), 
    percentage of frames, and lambda value.'''
    system = mda.Universe(f'bound/{run_folder}/vanish/input/SYSTEM.pdb', f'bound/{run_folder}/vanish/output/lambda-{lam_val}/traj000000001.dcd')
    no_close_waters = []
    no_frames_to_use = int(round(len(system.trajectory)*(percent_frames/100),0))
    first_frame = len(system.trajectory) - no_frames_to_use
    for frame in range(first_frame, len(system.trajectory)):
        system.trajectory[frame]
        no_close_waters.append(len(system.select_atoms(f'resname WAT and sphzone {length} index {index}'))/3)
    avg_close_waters = sum(no_close_waters)/len(no_close_waters)
    return round(avg_close_waters, 2)

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import os
    import sys
    import MDAnalysis as mda
    import numpy as np
    import re

    percent = 83
    index = 5192
    distance = 8 #A

    runs = int(sys.argv[1]) # only works for <10 runs

    run_folders = []

    # Get run dirs and ignore comments
    for run_no in range(1,runs+1):
        dir_names = os.listdir("bound")
        r = re.compile(f"run00{run_no}")
        dir_name = list(filter(r.match, dir_names))[0]
        run_folders.append(dir_name)

    # Define lambda values
    lam_folders = [x for x in os.listdir("bound/run001/vanish/output") if 'lambda' in x]
    lam_folders.sort()
    lam_vals = [x[-5:] for x in lam_folders]

    # Calculate average number of waters within 5 A of selected residue
    # of time percentages
    avg_close_waters_all = []
    for run_folder in run_folders:
        avg_close_waters = []
        for lam_val in lam_vals:
            avg_close_waters.append(avg_water(run_folder,lam_val,percent,index,distance))
        avg_close_waters_all.append(avg_close_waters)

    # Plot
    fig, ax = plt.subplots()
    for run_idx in range(len(run_folders)):
        ax.plot(lam_vals, avg_close_waters_all[run_idx], label=f"Run00{run_idx+1}")
    ax.set_xlabel(r'$\lambda$')
    ax.set_ylabel(f'Average No Waters\n within {distance} $\AA$ of index {index}')
    ax.legend()
    rounded_lam_vals =  [round(float(x),2) for x in lam_vals]

    ax.set_xticklabels(rounded_lam_vals)
    every_nth = 4
    for n, label in enumerate(ax.xaxis.get_ticklabels()):
        if n % every_nth != 0:
            label.set_visible(False)

    fig.tight_layout()
    fig.savefig(f'avg_waters_within_{distance}_index_{index}_percent_{percent}.png', dpi=1000, facecolor='w')
