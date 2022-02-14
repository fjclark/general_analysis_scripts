# Functions to get histograms of Boresch degrees of freedom over
# trajectories for specified value of lambda

import MDAnalysis as mda
from MDAnalysis.analysis.distances import dist
from MDAnalysis.lib.distances import calc_dihedrals
import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
from ..get_data import dir_paths
from ast import literal_eval
from ..save_data import mkdir_if_required


def get_mda_universe(leg, run, stage, lam_val):
    """Generate mda universe object.

    Args:
        leg (str): Bound or free
        run (i): Run number
        stage (str): Restrain, discharge, vanish
        lam_val (str): e.g. "0.000"

    Returns:
        mda.universe: mda universe object
    """
    run_name = dir_paths.get_run_name(run,leg)
    paths = dir_paths.get_dir_paths([run],leg)
    top_file = f"{paths[run_name][stage]['input']}/SYSTEM.top"
    traj_file = f"{paths[run_name][stage]['output']}/lambda-{lam_val}/traj000000001.dcd"
    u = mda.Universe(top_file, traj_file)
    print(f"Opening {traj_file}")

    return u

def read_boresch_rest(cfg_file):
    """Read config file to extract Boresch restraints dict.

    Args:
        cfg_file (str): Path to config file

    Returns:
        dict: Boresch restraints dict
    """
    with open(cfg_file,"r") as istream:
        lines = istream.readlines()

    boresch_dict = {}
    for l in lines:
        if l.startswith("boresch restraints dictionary"):
            dict_as_list = l.split("=")[1][1:-1] # remove leading space and \n
            boresch_dict = literal_eval("".join(dict_as_list))
            break

    return boresch_dict


# Functions to get Boresch DOF

def get_distance(idx1, idx2, u):
    """Distance between two atoms in Angstrom

    Args:
        idx1 (int): Index of first atom
        idx2 (int): Index of second atom
        u (mda universe): System

    Returns:
        float: Distance in Angstrom
    """
    distance = dist(mda.AtomGroup([u.atoms[idx1]]), mda.AtomGroup([u.atoms[idx2]]), box=u.dimensions)[2][0]
    return distance

def get_angle(idx1, idx2, idx3, u):
    """Angle between three particles in rad.

    Args:
        idx1 (int): Index of first atom
        idx2 (int): Index of second atom
        idx3 (int): Index of third atom
        u (mda universe): System

    Returns:
        float: Angle in rad
    """
    C = u.atoms[idx1].position 
    B = u.atoms[idx2].position 
    A = u.atoms[idx3].position 
    BA = A - B
    BC = C - B
    angle = np.arccos(np.dot(BA, BC)/(norm(BA)*norm(BC)))
    return angle

def get_dihedral(idx1, idx2, idx3, idx4, u):
    """Get dihedral based on four atom positions

    Args:
        idx1 (int): Index of first atom
        idx2 (int): Index of second atom
        idx3 (int): Index of third atom
        idx4 (int): Index of fourth atom
        u (mda universe): System

    Returns:
        float: Dihedral angle in rad
    """
    positions =[u.atoms[idx].position for idx in [idx1,idx2,idx3,idx4]]
    dihedral = calc_dihedrals(positions[0], positions[1], positions[2], positions[3], box = u.dimensions)
    return dihedral

def get_boresch_dof(anchor_ats, u):
    """Calculate the degrees of freedom defined by the Boresch
    restraints. In addition, get the internal bond angles made 
    by the anchor points in the receptor and ligand, because these
    are important for the stability of the dihedrals. Ordering of
    connection of anchor atoms is r3, r2, r1, l1, l2, l3.

    Args:
        anchor_ats (tuple): Tuple of anchor atom indices of 
        form (r1,r2,r3,l1,l2,l3), where r are anchors on the receptor
        and l on the ligand
        u (mda universe): System

    Returns:
        int, floats: Boresch degrees of freedom
    """
    r1,r2,r3,l1,l2,l3 = anchor_ats
    r = get_distance(r1,l1,u)
    thetaA = get_angle(r2,r1,l1,u)
    thetaB = get_angle(r1,l1,l2,u)
    phiA = get_dihedral(r3,r2,r1,l1,u)
    phiB = get_dihedral(r2,r1,l1,l2,u)
    phiC = get_dihedral(r1,l1,l2,l3,u)
    # Not restrained but distance from coolinearity must be checked
    thetaR = get_angle(r3,r2,r1,u) # Receptor internal angle
    thetaL = get_angle(l1,l2,l3,u) # Ligand internal angle
    return r, thetaA, thetaB, phiA, phiB, phiC, thetaR, thetaL


def track_dof(anchor_ats, u, percent_traj):
    """Get values, mean, and standard deviation of Boresch
    degrees of freedom and internal angles defined by supplied
    anchor atoms. Also calculate total variance accross all DOF
    , neglecting the internal angles.

    Args:
        anchor_ats (tuple): Anchor atom indices, of form (r1,r2,r3,l1,l2,l3)
        u (mda universe): The system
        percent_traj (float): Percentage of run to average over (25 % would
        result in intial 75 % of trajectory being discarded)

    Returns:
        dict: dictionary of form {"tot_var": tot_var, dof1 :{"mean":mean, "sd":sd,
         "values":[...]}, dof2:{...},...}
    """
    
    r1,r2,r3,l1,l2,l3 = anchor_ats
    n_frames = len(u.trajectory)
    first_frame = round(n_frames - ((percent_traj/100)*n_frames)) # Index of first frame to be used for analysis
    
    dof_dict = {}
    dof_list = ["r","thetaA","thetaB","phiA","phiB","phiC","thetaR","thetaL"]
    # Add sub dictionaries for each Boresch degree of freedom
    for dof in dof_list:
        dof_dict[dof]={}
        dof_dict[dof]["values"]=[]

    for i, frame in enumerate(u.trajectory):
        if i >= first_frame:
            if i == first_frame:
                print(f"First frame no: {i+1}")
            r, thetaA, thetaB, phiA, phiB, phiC, thetaR, thetaL = get_boresch_dof(anchor_ats,u)
            dof_dict["r"]["values"].append(r)
            dof_dict["thetaA"]["values"].append(thetaA)
            dof_dict["thetaB"]["values"].append(thetaB)
            dof_dict["phiA"]["values"].append(phiA)
            dof_dict["phiB"]["values"].append(phiB)
            dof_dict["phiC"]["values"].append(phiC)
            dof_dict["thetaR"]["values"].append(thetaR)
            dof_dict["thetaL"]["values"].append(thetaL)

            if i == n_frames-1:
                print(f"Last frame no: {i+1}")
                dof_dict["tot_var"]=0
                for dof in dof_list:
                    dof_dict[dof]["values"]=np.array(dof_dict[dof]["values"])
                    dof_dict[dof]["mean"]=dof_dict[dof]["values"].mean()
                    # For dihedrals, compute variance based on list of values 
                    # corrected for periodic boundary at pi radians, because there
                    #  is no problem with dihedrals in this region
                    if dof[:3] == "phi":
                        mean = dof_dict[dof]["mean"]

                        # correct variance - fully rigorous
                        corrected_values_sd = []
                        for val in dof_dict[dof]["values"]:
                            dtheta = abs(val - mean)
                            corrected_values_sd.append(min(dtheta, 2*np.pi-dtheta))
                        corrected_values_sd = np.array(corrected_values_sd) 
                        dof_dict[dof]["sd"]=corrected_values_sd.std()

                        # Correct mean (not exact and will fail if very well split above 
                        # and below 2pi)get middle of interval based on current mean
                        print("WARNING: Mean for dihedrals may be incorrect if true mean" \
                              " is near periodic boundary")
                        corrected_values_mean=[]
                        periodic_bound = mean - np.pi
                        if periodic_bound < -np.pi:
                            periodic_bound+=2*np.pi
                        # shift vals from below periodic bound to above
                        for val in dof_dict[dof]["values"]:
                            if val < periodic_bound:
                                corrected_values_mean.append(val+2*np.pi)
                            else:
                                corrected_values_mean.append(val)
                        corrected_values_mean = np.array(corrected_values_mean)
                        mean_corrected = corrected_values_mean.mean()
                        #shift mean back to normal range
                        if mean_corrected > np.pi:
                            dof_dict[dof]["mean"]=mean_corrected-2*np.pi
                        else:
                            dof_dict[dof]["mean"]=mean_corrected
                            
                    else:
                        dof_dict[dof]["sd"]=dof_dict[dof]["values"].std()
                    # Exclude variance of internal angles as these are not restrained
                    if (dof != "thetaR" and dof != "thetaL"):
                        dof_dict["tot_var"]+=dof_dict[dof]["sd"]**2
                    # Assume Gaussian distributions and calculate "equivalent"
                    # force constants for harmonic potentials
                    # so as to reproduce these distributions
                    dof_dict[dof]["k_equiv"]=0.593/(dof_dict[dof]["sd"]**2) # RT at 298 K is 0.593 kcal mol-1
    
    return dof_dict


def get_dof_dicts(leg, runs, stage, lam_val, percent_traj):
    """Get dof_dicts for given stage and lambda value
    for all supplied runs

    Args:
        leg (str): bound
        runs (list): Run numbers (ints)
        stage ([type]): restrain, discharge, or vanish
        lam_val (str): Window of interest
        percent_traj (float): Percentage of run to average over (25 % would
        result in intial 75 % of trajectory being discarded)

    Returns:
        dict: dict of dof_dicts, with run names as keys
    """
    if type(lam_val == float):
        lam_val = f"{lam_val:.3f}" # Account for input of float instead of string
    dof_dicts = {}
    paths = dir_paths.get_dir_paths(runs, leg)

    for run in runs:
        run_name = dir_paths.get_run_name(run,leg)
        dof_dicts[run_name] = {}
        u = get_mda_universe(leg, run, stage, lam_val)
        cfg_path = f'{paths[run_name][stage]["input"]}/sim.cfg'
        boresch_dict = read_boresch_rest(cfg_path)
        anchor_ats = tuple([x for x in boresch_dict["anchor_points"].values()])
        dof_dict = track_dof(anchor_ats, u, percent_traj)
        dof_dicts[run_name]=dof_dict

    return dof_dicts
    

def plot_dof_hists(leg, runs, stage, lam_val, percent_traj, selected_dof_list):
    """Plot histograms of selected degrees of freedom over specified
    final percentage of trajectory for supplied runs and lambda window.

    Args:
        leg (str): bound
        runs (list): [
        runs (list): Run numbers (ints)
        stage ([type]): restrain, discharge, or vanish
        lam_val (str): Window of interest
        percent_traj (float): Percentage of run to average over (25 % would
        result in intial 75 % of trajectory being discarded)
        selected_dof_list (list): Subset of ["r","thetaA","thetaB","phiA","phiB",
        "phiC","thetaR","thetaL"]
    """

    dof_dicts = get_dof_dicts(leg, runs, stage, lam_val, percent_traj)
    no_dof = len(selected_dof_list)

    fig, axs = plt.subplots(1, no_dof, figsize=(4*no_dof,4), dpi=1000)
    colours =  ['#00429d', '#3a8bbb', '#ffb59a', '#ff6b95', '#93003a'] # Will cause issues for more than 5 runs

    for j, run in enumerate(runs):
        run_name = dir_paths.get_run_name(run,leg)
        for i, dof in enumerate(selected_dof_list):
            ax = axs[i]
            values = dof_dicts[run_name][dof]["values"]
            mean = dof_dicts[run_name][dof]["mean"]
            sd = dof_dicts[run_name][dof]["sd"]
            ax.hist(values,label = f"{run_name}", color=colours[j], edgecolor='k')
            ax.axvline(mean, linestyle = "dashed", color=colours[j], linewidth=2, label=f"Mean: {mean:.2f}\nSD: {sd:.2f}")
            if dof == "r":
                ax.set_xlabel("r ($\AA$)")
            else:
                ax.set_xlabel(f"{dof} (rad)")
            ax.set_ylabel("Counts")
            ax.legend(loc=(1.04,0))

    fig.tight_layout()
    mkdir_if_required("analysis/comparitive_analysis")
    fig.savefig(f"analysis/comparitive_analysis/{leg}_{stage}_{lam_val}_boresch_dof_hists.png")


def plot_dof_vals(leg, runs, stage, lam_val, percent_traj, selected_dof_list,legend):
    """Plot values of selected degrees of freedom over specified
    final percentage of trajectory for supplied runs and lambda window.

    Args:
        leg (str): bound
        runs (list): [
        runs (list): Run numbers (ints)
        stage ([type]): restrain, discharge, or vanish
        lam_val (str): Window of interest
        percent_traj (float): Percentage of run to average over (25 % would
        result in intial 75 % of trajectory being discarded)
        selected_dof_list (list): Subset of ["r","thetaA","thetaB","phiA","phiB",
        "phiC","thetaR","thetaL"]
    """

    dof_dicts = get_dof_dicts(leg, runs, stage, lam_val, percent_traj)
    no_dof = len(selected_dof_list)

    fig, axs = plt.subplots(1, no_dof, figsize=(4*no_dof,4), dpi=1000)
    colours =  ['#00429d', '#3a8bbb', '#ffb59a', '#ff6b95', '#93003a'] # Will cause issues for more than 5 runs

    for j, run in enumerate(runs):
        run_name = dir_paths.get_run_name(run,leg)
        for i, dof in enumerate(selected_dof_list):
            ax = axs[i]
            values = dof_dicts[run_name][dof]["values"]
            mean = dof_dicts[run_name][dof]["mean"]
            sd = dof_dicts[run_name][dof]["sd"]
            ax.plot([x for x in range(len(values))], values, label = f"{run_name}", color=colours[j])
            ax.axhline(mean, linestyle = "dashed", color=colours[j], linewidth=2, label=f"Mean: {mean:.2f}\nSD: {sd:.2f}")
            if dof == "r":
                ax.set_ylabel("r ($\AA$)")
            else:
                ax.set_ylabel(f"{dof} (rad)")
            ax.set_xlabel("Frame No")
            if legend:
                ax.legend(loc=(1.04,0))

    fig.tight_layout()
    mkdir_if_required("analysis/comparitive_analysis")
    fig.savefig(f"analysis/comparitive_analysis/{leg}_{stage}_{lam_val}_boresch_dof_vals.png")


# Just seems to return whitespace (but no errors) if functions modified to return figures
#def plot_dof(leg, runs, stage, lam_val, percent_traj, selected_dof_list):
#    fig = plt.figure(figsize=(4*len(selected_dof_list), 8))
#    subfigs = fig.subfigures(1, 2)
#    subfigs[0] = plot_dof_hists(leg, runs, stage, lam_val, percent_traj, selected_dof_list)
#    subfigs[1] = plot_dof_vals(leg, runs, stage, lam_val, percent_traj, selected_dof_list, False)
#    fig.savefig(f"analysis/{leg}_{stage}_{lam_val}_boresch_dof.png")