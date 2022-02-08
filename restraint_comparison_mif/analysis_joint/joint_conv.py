# Plot comparitive convergence (overall and for individual stages) for
# all runs in conv_dict, which is of the form {run:{stage:{cumtime1:
# {dg_tot:DG_tot, pmf:{lamval1: DG1, lamval2:DG2, ...}}, cumtime2:{...}, ...}}} 
# Assumes sample size of 5

import matplotlib.pyplot as plt
import pickle
import numpy as np

LEG = "bound"


def plot_conv(ax, leg, stage, x_vals, y_vals_list):
    """Plot convergence with cumulative sampling time 
    for a given stage.

    Args:
        ax (ax): Axis on which to plot
        leg (str): Bound or free
        stage (stage): Restrain, discharge, or vanish
        x_vals (list): Cumulative sampling times
        y_vals_list (list): List of lists of total dg with cumulative sampling time
    """
    y_vals = np.array(y_vals_list)
    y_avg = np.mean(y_vals, axis=0)
    conf_int = np.std(y_vals, axis=0)*(2.776/np.sqrt(5)) # 95% C.I. for sample size of 5, so t = 2.776 (4 DOF)

    ax.plot(x_vals,y_avg)
    for i, entry in enumerate(y_vals):
        ax.plot(x_vals,entry, linestyle='dashed', label=f"run {i+1}")
    ax.set_title(f'{leg} {stage}')
    ax.set_ylabel('$\Delta \it{G}$ / kcal.mol$^-$$^1$')
    ax.set_xlabel('Cumulative sampling time / ns')
    ax.legend()
    ax.fill_between(x_vals, y_avg-conf_int, y_avg+conf_int, alpha=0.5, facecolor='#ffa500')


def plot_stages_conv(pickled_data):
    """Plot convergence accross all runs for individual stages.

    Args:
        pickled_data (str): Path to pickled convergence dictionary data
    """
    with open(pickled_data, "rb") as istream:
        conv_dict = pickle.load(istream)

    runs = list(conv_dict.keys())
    stages = list(conv_dict[runs[0]].keys())
    no_stages = len(stages)

    fig, axs = plt.subplots(1, no_stages, figsize = (4*no_stages,4), dpi = 1000)
    for i, stage in enumerate(stages):
        cumtimes = list(conv_dict[runs[0]][stage].keys())
        fr_nrgs = []
        for run in runs:
            fr_nrg_run = []
            for cumtime in cumtimes:
                fr_nrg_run.append(conv_dict[run][stage][cumtime]["dg_tot"])
            fr_nrgs.append(fr_nrg_run)

        plot_conv(axs[i], LEG, stage, cumtimes, fr_nrgs)

    fig.tight_layout()
    fig.savefig(f"analysis/{LEG}_stages_joint_convergence.png")

# def plot_overall_conv(pickled_data)

if __name__ == "__main__":
    import sys
    data = sys.argv[1]
    plot_stages_conv(data)
    #plot_overall_conv(data)