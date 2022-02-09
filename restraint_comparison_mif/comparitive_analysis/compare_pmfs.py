# Plot all PMFs for a given leg (at max cumulative sampling time) on the same axes.

from ..get_data.convergence_data import read_mbar_data
from ..get_data.dir_paths import get_dir_paths
from .compare_conv import plot_conv
import matplotlib.pyplot as plt
import numpy as np

def plot_all_pmfs(run_nos, leg):
    """Plot PMFs for all runs for individual stages. Includes 95% CI.

    Args:
        run_nos (list): List of run numbers to plot
        leg (str): Bound or free
    """
    dir_paths = get_dir_paths(run_nos, leg)
    runs = list(dir_paths.keys())
    stages = list(dir_paths[runs[0]].keys()) # stages from first run
    no_stages = len(stages)

    fig, axs = plt.subplots(1, no_stages, figsize = (4*no_stages,4), dpi = 1000)
    for i, stage in enumerate(stages):
        lam_vals = dir_paths[runs[0]][stage]["lam_vals"]
        pmfs = []
        for run in runs:
            mbar_path = dir_paths[run][stage]["mbar_data"]
            _, pmf, _ = read_mbar_data(mbar_path, lam_vals) # Throw away overall DG and overlap
            pmfs.append(pmf)

        plot_conv(axs[i], leg, stage, lam_vals, pmfs, "$\lambda$", "$\Delta \it{G}$ / kcal.mol$^-$$^1$")

    fig.tight_layout()
    fig.savefig(f"analysis/{leg}_pmf_comparison")


if __name__ == "__main__":
    import sys
    run_nos = sys.argv[1]
    leg = sys.argv[2]
    plot_all_pmfs(run_nos, leg)