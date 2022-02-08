# Based on convergence-plots.py (see https://github.com/michellab/MDM2-DG_paper)
# This produces convergence data for the stages specified (both overall and for
# individual stages). The resulting dictionary is of the form:
# {run:{stage:{cumtime1: {dg_tot:DG_tot, pmf:{lamval1: DG1, lamval2:DG2, ...}}, cumtime2:{...}, ...}}}

from . import get_dir_paths
import os, pickle

#import glob,sys,os
#import numpy as np
#import math
#import matplotlib.pyplot as plt

# Constants TODO: Add arg parser
CHUNKSIZE = 1  # ns - how far to separate each calculation of energies
TIMESTEP = 0.000004  # ns
NRGFREQ = 100  # How many steps between energy evaluations
LEG = "bound"
RUN_NOS = [1, 2, 3, 4, 5]
SIMTIME = {"restrain": {"wind_len": 6, "discard": 1}, "discharge": {
    "wind_len": 6, "discard": 1}, "vanish": {"wind_len": 8, "discard": 3}}
DIR_PATHS = get_dir_paths.get_dir_paths(RUN_NOS, LEG)


def truncate_simfile(in_file, out_file, start_time, end_time):
    """Truncate a simfile between given start and end time.

    Args:
        in_file (str): Path to input file
        out_file (str): Path to input file 
        start_time (int): Start time in ns
        end_time (int): End time in ns
    """
    with open(in_file, "r") as istream:
        with open(out_file, "w") as ostream:
            for line in istream.readlines():
                if line.startswith("#"):
                    ostream.write(line)
                    continue
                elems = line.split()
                time = (float(elems[0])+NRGFREQ)*TIMESTEP
                if (time < start_time):
                    continue
                if (time > end_time):
                    break
                ostream.write(line)


def do_mbar(lam_vals, input_dir, output_dir, start_time, end_time):
    """Perform MBAR analysis for supplied lambda values for data between
    given start and end times. 

    Args:
        lam_vals (list): List of lambda values for which MBAR will be performed
        input_dir (str): Output directory for desired stage
        output_dir (str): Directory where output will be saved
        start_time (float): Start time in ns
        end_time (float): End time in ns

    Returns:
        (float): cumulative sampling time 
    """
    cumtime = 0
    delta_t = end_time - start_time

    for lam_val in lam_vals:
        os.system(f"mkdir {output_dir}/lambda-{lam_val}")
        truncate_simfile(f"{input_dir}/lambda-{lam_val}/simfile.dat",
                         f"{output_dir}/lambda-{lam_val}/simfile.dat", start_time, end_time)
        cumtime += delta_t

    cmd = f"/home/finlayclark/anaconda3/envs/biosimspace-dev/bin/analyse_freenrg mbar \
         -i {output_dir}/lambda*/simfile.dat -p 100 --temperature 298.0 > {output_dir}/mbar.dat"
    os.system(cmd)

    return cumtime


def read_mbar_data(mbar_file, lam_vals):
    """Read the total MBAR free energy estimate and MBAR
    PMF values

    Args:
        mbar_file (str): Path to .dat file to read
        lam_vals (list): List of lam vals in string format

    Returns:
        dg_tot (float): Total MBAR free energy for stage
        pmf_dict (dict): Dict of PMF values from MBAR for each lambda value
    """
    dg_tot = 0
    pmf_dict = {}

    with open(mbar_file,"r") as istream:
        lines = istream.readlines()

        pmf_idx = 0
        dg_tot_idx = 0
        for i, line in enumerate(lines):
            if line.startswith("#PMF from MBAR"):
                pmf_idx = i+1
            if line.startswith("#MBAR free energy difference"):
                dg_tot_idx = i+1

        for i, lam_val in enumerate(lam_vals):
            pmf_dict[lam_val] = float(lines[pmf_idx+i].split(" ")[1].strip())
        dg_tot = float(lines[dg_tot_idx].split()[0].strip(","))
            
    return dg_tot, pmf_dict


def get_convergence_dict():
    """Create dictionary of the form {run:{stage:{cumtime1: {dg_tot:DG_tot, 
    pmf:{lamval1: DG1, lamval2:DG2, ...}}, cumtime2:{...}, ...}}} to store all
    convergence data.

    Returns:
        dict: convergence data
    """

    conv_dict = {}
    
    for run in DIR_PATHS.keys():
        conv_dict[run]={}
        for stage in DIR_PATHS[run].keys():
            conv_dict[run][stage]={}
            start_time = SIMTIME[stage]["discard"]
            final_end_time = SIMTIME[stage]["wind_len"]
            end_times = [(x*CHUNKSIZE)+start_time for x in range(1, int((final_end_time-start_time)/CHUNKSIZE)+1)]
            cum_times = [x-start_time for x in end_times]

            for i, cumtime in enumerate(cum_times):
                conv_dict[run][stage][cumtime]={}
                lam_vals = DIR_PATHS[run][stage]["lam_vals"]
                input_dir = DIR_PATHS[run][stage]["output"] # output of simulations is input to do_mbar()

                os.system("mkdir tmp")
                print("#############################################################################################")
                print(f"MBAR analysis in progress for {run}, {stage}, cumulative sampling time {cumtime} ns")
                print("#############################################################################################")
                do_mbar(lam_vals, input_dir, "./tmp", start_time,end_times[i])
                dg_tot, pmf_dict = read_mbar_data("./tmp/mbar.dat",lam_vals)
                conv_dict[run][stage][cumtime]["dg_tot"] = dg_tot
                conv_dict[run][stage][cumtime]["pmf"] = pmf_dict
                print(f"dg_tot = {dg_tot}")
                print(f"pmf_dict = {pmf_dict}")
                os.system("rm -rf tmp")

    if not os.path.isdir("./analysis"):
        os.system("mkdir analysis")
    dumpme = open("analysis/convergence_data.pickle","wb") # save data to file
    pickle.dump(conv_dict, dumpme)
    dumpme.close()

    return conv_dict

if __name__ == "__main__":
    print(get_convergence_dict())