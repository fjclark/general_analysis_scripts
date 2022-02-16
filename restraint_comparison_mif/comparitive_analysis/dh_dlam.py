# Plot the gradient dH/dlam with time for all windows

import matplotlib.pyplot as plt
from ..get_data.dir_paths import get_dir_paths
from ..save_data import mkdir_if_required


def read_simfile(simfile_path, percent_traj, timestep=4, nrg_freq =100):
    """Read gradients over specified percentage of end of trajectory from
    simfile.

    Args:
        simfile_path (str): Path to simfile.dat
        percent_traj (float): Percentage of data to use, beginning from end of simulation
        timestep (int, optional): Timestep in fs. Defaults to 4 fs.
        nrg_freq (int, optional): Frequency of energy saving. Required because this changes 
        form of the simfile. Defaults to 100.

    Returns:
        times (list): The times from the specified percentage of the trajectory
        grads (list): The gradients at the given times.
    """
    steps = []
    grads = []

    with open(simfile_path, "rt") as f:
        lines = f.readlines()

    for l in lines:
        vals = l.split()
        if not l.startswith("#"):
            step = int(vals[0].strip())
            grad = float(vals[2].strip())
            steps.append(step)
            grads.append(grad)

    times = [(x+nrg_freq)*timestep/1000000 for x in steps] # Convert steps to time in ns
    n_steps = len(steps)
    first_step_idx = round(n_steps - ((percent_traj/100)*n_steps)) # Index of first step to be used for analysis
    print(f"Reading simfile from {times[first_step_idx]} ns to {times[-1]} ns")

    times = times[first_step_idx:]
    grads = grads[first_step_idx:]

    return times, grads


def plot_grad(run_name, simfile_path, ax, percent_traj, timestep =4, nrg_freq = 100):
    """Plot gradients with time on supplied axis for given run, leg, stage,
    and lambda value.

    Args:
        simfile_path (str): Path to simfile
        ax (matplotlib axis): axis on which to plot 
        percent_traj (float): Percentage of end of trajectory over which to 
        timestep (int, optional): Timestep in fs. Defaults to 4 fs.
        nrg_freq (int, optional): Frequency of energy saving. Required because this changes 
        form of the simfile. Defaults to 100.
    """
    times, grads = read_simfile(simfile_path, percent_traj, timestep, nrg_freq)
    
    ax.plot(times, grads, label=f"Run {run_name}", alpha=0.3)
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("d$H$/d$\lambda$ / kcal.mol$^-$$^1$")


def plot_grads(leg, runs, percent_traj_dict, timestep=4, nrg_freq=100):
    """Plot gradients against time for all supplied runs, for 
    all lambda windows in all stages of given leg.

    Args:
        run_name (str): Name of run
        leg (str): bound or free
        runs (list): List of run numbers (ints)
        percent_traj_dict (dict): Specifies percentage of trajectories
        to use for given stage (of form {"stage":percent,...})
        timestep (int, optional): Timestep in fs. Defaults to 4 fs.
        nrg_freq (int, optional): Frequency of energy saving. Required because this changes 
        form of the simfile. Defaults to 100.
    """
    paths = get_dir_paths(runs, leg)
    run_names = list(paths.keys())
    stages = list(paths[run_names[0]].keys())
    stages = ["vanish"] + stages # split vanish into two - easiest to add another vanish
    fig, _ = plt.subplots(1,1,figsize=(20,52), dpi=200)
    subfigs = fig.subfigures(1,1*len(stages)) # subfigs to allow different no plots in each column

    vanish_count = 0 # split vanish into two stages and count number of vanish stages already plotted

    for i, stage in enumerate(stages):
        lam_vals = paths[run_names[0]][stage]["lam_vals"]
        if stage == "vanish": # split in half as many windows
            median_idx = round(len(lam_vals)/2)
            if vanish_count == 0:
                lam_vals =lam_vals[:median_idx]
                vanish_count +=1
            else:
                lam_vals = lam_vals[median_idx:]

        subfig = subfigs[i]
        subfig.suptitle(f'{leg} {stage}')
        subfig_axs = subfig.subplots(len(lam_vals),1, sharex=True)
        percent_traj = percent_traj_dict[stage]

        for j, lam_val in enumerate(lam_vals):
            ax = subfig_axs[j]
            ax.title.set_text(f"$\lambda$ = {lam_val}")
            for run_name in run_names:
                print(f"Plotting grads for {leg}, {stage}, lambda: {lam_val}, {run_name}")
                lam_path = paths[run_name][stage]["lam_paths"][j]
                simfile_path = f"{lam_path}/simfile.dat"
                plot_grad(run_name, simfile_path, ax, percent_traj, timestep, nrg_freq)
            ax.legend(loc="best")
            # Keep only the bottom time 
            if j!= len(lam_vals) - 1:
                x_axis = ax.get_xaxis()
                x_axis.set_visible(False)
        
    mkdir_if_required("analysis/overall_convergence")
    fig.savefig(f'analysis/overall_convergence/{leg}_grads.png')