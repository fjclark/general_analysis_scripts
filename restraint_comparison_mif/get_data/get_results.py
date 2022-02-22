# Extracts all contributions to the free energy of binding for given leg (from MBAR)
# along with associated error calculated by propagation of the MBAR errors 
# (a large underestimation of the uncertainty calculated from multiple runs).

from cProfile import run
import numpy as np
import statsmodels.stats.api as sms

from .dir_paths import get_dir_paths
from ..save_data import mkdir_if_required
from .convergence_data import read_mbar_data
import scipy.stats as st


def get_lj_corr(input_file):
    """Finds LJ correction from .dat file

    Args:
        input_file (str): Input file

    Returns:
        : tuple of (LJ correction, standard deviation)
    """
    with open(input_file,'r') as file:
        lines = file.readlines()
    
    try:
        correction = float(lines[0].split()[2])
        conf_int= float(lines[0].split()[4])
    except:
        print(f"ERROR: Unable to read {input_file}. LJ correction has likely failed")
        correction = 0
        conf_int = 0
        
    return correction,conf_int


def get_boresch_ana(input_file):
    """Finds analytical correction for releasing Boresch 
    restraints and returns correction as string

    Args:
        input_file (str): Input file

    Returns:
        energy: analytical correction for releasing restraints
    """
    with open(input_file,'r') as file:
        lines = file.readlines()
        
    for i, line in enumerate(lines):
        if line.startswith('Analytical correction for releasing Boresch restraints'):
            energy=float(lines[i].split()[-3])

    conf_int = 0 # Include this for consistency with other functions, which return non-zero SD
            
    return energy, conf_int


def get_results(leg = "bound", run_nos = [1,2,3,4,5]):
    """Get the MBAR free energy estimates and standard deviations
    associated with all stages in a given leg for the supplied runs

    Args:
        leg (str, optional): Free or bound. Defaults to "bound".
        run_nos (list, optional): List of int run numbers. Defaults to [1,2,3,4,5].

    Returns:
        dict: Dictionary of form {"run001":{"vanish":dg, "restrain":dg, "lj_corr":dg,...} ,...}
    """

    paths = get_dir_paths(run_nos, leg)
    results = {}
    
    for run_name in paths.keys():
        results[run_name]={}
        for stage in paths[run_name].keys():
            mbar_file = paths[run_name][stage]["mbar_data"]
            lam_vals = paths[run_name][stage]["lam_vals"]
            dg, conf_int, _, _ = read_mbar_data(mbar_file,lam_vals) # throw away PMF and overlap
            results[run_name][stage] = (dg, conf_int)

            if stage == "vanish":
                output_dir = paths[run_name][stage]["output"]
                dg, conf_int = get_lj_corr(f"{output_dir}/freenrg-LJCOR.dat")
                results[run_name]["lj_corr"] = (dg, conf_int)
                # TODO: Modify this to work if not Boresch 
                if leg == "bound":
                    dg, conf_int = get_boresch_ana(f"{output_dir}/boresch_analytical_correction.dat")
                    results[run_name]["boresch_ana_corr"] = (dg, conf_int)
                    # Symmetry corrections assume 298 K (RT = 0.592187)
                    results[run_name]["symm_corr_binding_sites_298"] = (0.65, 0) # Three-fold symmetry of binding sites (so RTln3)
                    results[run_name]["symm_corr_phenol_298"] = (0.41, 0) # Rotation of phenol hindered in binding site (so RTln2)
            
        dg_tot = sum([x[0] for x in results[run_name].values()]) # Energy is first value in tuple
        var_tot = sum([x[1]**2 for x in results[run_name].values()])
        ci_tot = np.sqrt(var_tot)
        results[run_name]["dg_tot"] = (dg_tot, ci_tot) 
                    
    return results


def write_results_indiv(results):
    """Write results for individual runs to text file.

    Args:
        results (dict): Dictionary of results (possibly for several runs)
    """
    for run_name in results.keys():
        mkdir_if_required("analysis")
        mkdir_if_required("analysis/results")
        mkdir_if_required(f"analysis/results/{run_name}")

        with open(f"analysis/results/{run_name}/{run_name}_results.txt", "wt") as f:
            f.write("Free energy estimates from MBAR:\n")
            f.write("(Uncertainties are 95 % C.I.s derived from MBAR for single run)\n\n")
            for contribution in results[run_name].keys():
                dg = results[run_name][contribution][0]
                sd = results[run_name][contribution][1]
                f.write(f"{contribution}: {dg:.2f} +/- {sd:.2f} kcal/mol\n")


def write_results_overall(results):
    """Write summary results to text file. Calculate
    95 % C.I. for all stages.

    Args:
        results (dict): Results dictionary
    """
    mkdir_if_required("analysis")
    mkdir_if_required("analysis/results")

    tot_dict = {}
    run_names = list(results.keys())
    for contribution in results[run_names[0]]: # Use first run to check contributions to free energy
        tot_dict[contribution]={}
        tot_dict[contribution]["values"]=[]
        for run_name in run_names:
            tot_dict[contribution]["values"].append(results[run_name][contribution][0]) # Ignore SD from individual runs
        tot_dict[contribution]["values"] =np.array(tot_dict[contribution]["values"])

        vals = tot_dict[contribution]["values"]
        tot_dict[contribution]["dg"] = vals.mean()
        conf_int = st.t.interval(0.95, len(vals)-1, loc=np.mean(vals), scale=st.sem(vals))
        tot_dict[contribution]["95% C.I."] = vals.mean() - conf_int[0] # Because C.I. returned as tuple (min, max) 

    with open(f"analysis/results/results.txt", "wt") as f:
        f.write("Overall free energy estimates from MBAR:\n")
        f.write(f"(Uncertainty estimates are 95 % C.I.s derived from differences between {len(run_names)} replicate runs)\n\n")
        for contribution in tot_dict.keys():
            dg = tot_dict[contribution]["dg"]
            conf = tot_dict[contribution]["95% C.I."]
            f.write(f"{contribution}: {dg:.2f} +/- {conf:.2f} kcal/mol\n")


def write_results(leg="bound", run_nos = [1,2,3,4,5]):
    """Retrieve results for given leg for all supplied runs.
    Save individual summary and overall summary with 95 % C.I.s

    Args:
        leg (str, optional): Bound or free. Defaults to "bound".
        run_nos (list, optional): List of run numbers (ints). Defaults to [1,2,3,4,5].
    """
    print("###############################################################################################")
    print("Writing results")
    results = get_results(leg, run_nos)
    write_results_indiv(results)
    write_results_overall(results)