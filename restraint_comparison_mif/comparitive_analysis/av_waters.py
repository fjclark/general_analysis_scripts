# Scripts to plot average number of waters within certain distance
# of given reference for given stage (most likely bound vanish)

import matplotlib.pyplot as plt
from .boresch_dof import get_mda_universe
from .compare_conv import plot_conv
from ..get_data.dir_paths import get_dir_paths
from ..save_data import mkdir_if_required

def get_av_waters(u, percent_traj, index, length):
    """Calculate average number of waters within given distance
    of atom with given index over the specified percentage of the end
    of the trajectory.

    Args:
        u (mda universe): system and trajectory
        percent_traj (float): percentage of trajectory (beginning from end) over 
        which to average
        index (int): Atom fromm which distance is calculated
        length (float): Distance in Angstrom

    Returns:
        float: average number of waters
    """
    no_close_waters = []
    n_frames = len(u.trajectory)
    first_frame = round(n_frames - ((percent_traj/100)*n_frames)) # Index of first frame to be used for analysis
    print(f"Calculating av. no. waters within {length} A of index {index} ")

    for frame in range(first_frame, n_frames):
        u.trajectory[frame]
        no_close_waters.append(len(u.select_atoms(f'resname WAT and sphzone {length} index {index}'))/3) # /3 as returns all atoms in water
    avg_close_waters = sum(no_close_waters)/len(no_close_waters)

    return avg_close_waters


def plot_av_waters(leg="bound", runs=[1,2,3,4,5], stage="vanish", percent_traj=62.5, index=20, length=8):
    """Plot the average number of waters within the specified distance of
    a given atom over all lambda windows for all runs, using the specified 
    percentage of all trajectories

    Args:
        leg (str): bound or free
        runs (list): list of runs (ints)
        stage (str): e.g. "vanish"
        percent_traj (float): percentage of trajectory (beginning from end) over 
        which to average
        index (int): Atom fromm which distance is calculated
        length (float): Distance in Angstrom
    """

    paths = get_dir_paths(runs, leg)
    run_names = list(paths.keys())
    lam_vals = paths[run_names[0]][stage]["lam_vals"]

    av_waters_all_runs = [] # List of lists - av water at each lam val for each run

    for run_no in runs:
        av_waters =[]
        for lam_val in lam_vals:
            u = get_mda_universe(leg, run_no, stage, lam_val)
            waters = get_av_waters(u, percent_traj, index,length)
            av_waters.append(waters)
        av_waters_all_runs.append(av_waters)

    fig, ax = plt.subplots(figsize = (4,4), dpi=1000)
    plot_conv(ax, leg, stage, lam_vals, av_waters_all_runs, "$\lambda$",
              f"Av. no waters within \n{length} $\AA$ of index {index}")

    fig.tight_layout()
    mkdir_if_required("analysis/comparitive_analysis")
    fig.savefig(f'analysis/comparitive_analysis/{leg}_{stage}_{index}_{length}_av_waters.png')