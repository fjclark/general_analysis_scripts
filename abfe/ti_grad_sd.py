#!/bin/python

def get_grad_std(datafiles):
    """Returns array of average gradients from TI from list of datafiles.

    Args:
        datafiles (list): List of datafiles

    Returns:
        lam_vals (list): List of lambda values 
        ti_gradients_std (np.array): 1D array of std of average gradients for TI
        at each lam val
    """
    
    ti_gradients_list = []
    lam_vals = []

    # Get no lambda values and lambda values
    # Read the first data file
    with open(datafiles[0],'r') as f:
        lines = f.readlines()

    line_idx_ti_avg_grad = 0
    line_idx_pmf_ti = 0

    for i,l in enumerate(lines):
        if l.startswith('#TI average gradients and standard deviation in kcal/mol'):
            line_idx_ti_avg_grad = i
        if l.startswith('#PMF from TI in kcal/mol'):
            line_idx_pmf_ti = i

    for i in range(line_idx_ti_avg_grad+1, line_idx_pmf_ti):
        lam_vals.append(np.array(lines[i].strip().split(' '))[0].astype('float'))

    # Read through all files and add each set of gradients to a row in 
    for file in datafiles:
        with open(file, 'r') as f:
            lines=f.readlines()

        gradients = []

        for i in range(line_idx_ti_avg_grad+1, line_idx_pmf_ti):
            gradients.append(np.array(lines[i].strip().split(' '))[1].astype('float'))
            
        ti_gradients_list.append(gradients)

    # Array of gradients
    ti_gradients_array = np.stack(ti_gradients_list)

    # Get SD
    ti_gradients_sd = ti_gradients_array.std(axis=0)

    # Print output
    print(f"Lambda values and associated standard deviation:")
    for i, j in zip(lam_vals, ti_gradients_sd):
        print(f"Lam val: {i:.3f}, SD: {j:.3f}")

    return lam_vals, ti_gradients_sd


if __name__ == "__main__":

    import numpy as np
    import matplotlib.pylab as plt
    import seaborn as sbn
    import sys
    sbn.set(font_scale=2,style="white")

    # Produces plot of std of gradient from TI for bound vanish stage for 
    # specified number of runs

    runs = int(sys.argv[1])
    #runs = 3

    # get files
    datafiles = []
    for i in range(runs):
        datafiles.append(f'bound/run00{i+1}/vanish/output/freenrg-MBAR-p-83-overlap.dat')

    # Get grad std
    lam_vals, ti_gradients_std = get_grad_std(datafiles)

    # Plot
    fig, ax = plt.subplots(figsize=(12,12))
    ax.set_xlabel('$\lambda$')
    ax.set_ylabel(f'SD of Avg Grad over {runs} Runs [kcal / mol]')
    ax.plot(lam_vals, ti_gradients_std)
    for index in range(len(lam_vals)):
        ax.text(lam_vals[index],ti_gradients_std[index],lam_vals[index],size=12)
    fig.tight_layout()
    fig.savefig(f'ti_sd_{runs}_runs.png')