# Run complete set of default analyses 
import os

from .get_data import convergence_data
from .get_data import get_results
from .comparitive_analysis import compare_conv
from .comparitive_analysis import compare_pmfs
from .comparitive_analysis import overlap
from .comparitive_analysis import indiv_pmf_conv
from .comparitive_analysis import dh_dlam
from .comparitive_analysis import restrained_dof
from .comparitive_analysis import rmsd
from .comparitive_analysis import av_waters

def run_analysis(leg = "bound", run_nos=[1,2,3,4,5], restraint_type="Boresch", timestep=4, nrg_freq=100,
percent_traj_dict = {"restrain":83.33333333, "discharge":83.33333333, "vanish":62.5, "rigidify":83.33333333,
                     "unrigidify_lig":83.33333333, "unrigidify_prot":83.33333333},
                      simtime = {"restrain": {"wind_len": 6, "discard": 1}, "discharge": {
                          "wind_len": 6, "discard": 1}, "vanish": {"wind_len": 8, "discard": 3},
                          "release": {"wind_len": 2, "discard": 1}, "unrigidify_lig": {
                          "wind_len": 6, "discard": 1},"unrigidify_prot": {
                          "wind_len": 6, "discard": 1}, "rigidify": {
                          "wind_len": 6, "discard": 1}, "release_2": {"wind_len": 2, "discard": 1}}
                    ):

                    # For free leg, remember to change to:
                    #percent_traj_dict = {"discharge":83.33333333, "vanish":83.3333333}):
                    #simtime = {"discharge": {"wind_len": 6, "discard": 1}, "vanish": {"wind_len": 6, "discard": 1}} # ns

    """Run analysis of bound leg. Note that some global variables must be changed in convergence_data.py (this will
    be fixed in future). If the convergence data has already been generated, the analysis will start from there.

    Args:
        leg (str, optional): Does not currently allow 'free' as an option. Defaults to "bound".
        run_nos (list, optional): _description_. Defaults to [1,2,3,4,5].
        restraint_type (str, optional): Boresch, multiple_dist, or Cart. Defaults to "Boresch".
        timestep (int, optional): In fs. Defaults to 4.
        nrg_freq (int, optional): Steps between energy evaluations. Defaults to 100.
        percent_traj_dict (dict, optional): Percentage of trajectory to use for analysis for each stage.
        Defaults to {"restrain":83.33333333, "discharge":83.33333333, "vanish":62.5}.
        simtime (dict): Lengths of simulations and lengths of inital periods to discard as equilibration, in ns
                        , for generation of the convergence data. 
    """

    print("###############################################################################################")
    print(f"Analysing the {leg} leg for runs: {run_nos} and calculation type = {restraint_type}")
    print("Ensure you are in the base directory and the development version of biosimspace is activated")

    # Only calculate convergence data if this has not been done already
    if not os.path.isfile("analysis/convergence_data.pickle"):
        convergence_data.get_convergence_dict(leg=leg, run_nos=run_nos, nrg_freq=nrg_freq,
                                              timestep=timestep/1000000, # Convert to ns
                                               simtime=simtime)

    if leg == "bound":
        get_results.write_results(leg, run_nos, restraint_type)
        compare_conv.plot_stages_conv("analysis/convergence_data.pickle", leg)
        compare_conv.plot_overall_conv("analysis/convergence_data.pickle", leg)
        compare_pmfs.plot_all_pmfs(run_nos, leg)
        overlap.plot_overlap_mats(leg, run_nos)
        indiv_pmf_conv.plot_pmfs_conv(leg, run_nos)
        dh_dlam.plot_grads(leg, run_nos, percent_traj_dict, timestep, nrg_freq)

        # Plot average waters within 8 A of CG2 in VAL and N in PRT on opposite sides of binding pocket.
        # This gives reasonable coverage of the pocket while excluding most waters outside.
        av_waters.plot_av_waters(leg, run_nos, stage="vanish", 
        percent_traj=percent_traj_dict["vanish"], index=1637,length=8, index2=34, length2=8)
        # Plot DOF for restrain lam = 0, restrain lam = 1, discharge lam = 1 and vanish lam =1
        plot_winds = [("restrain",0.000),("restrain",1.000),("discharge",1.000),("vanish",1.000)]

        if restraint_type == "Boresch":
            selected_dof_list = ["r","thetaA","thetaB","phiA","phiB","phiC"]
        elif restraint_type == "Cart":
            selected_dof_list = ["xr_l1", "yr_l1", "zr_l1", "phi", "theta", "psi"]
        elif restraint_type == "multiple_dist":
            selected_dof_list = []

        for wind in plot_winds:
            restrained_dof.plot_dof_hists(leg, run_nos, wind[0], wind[1], percent_traj_dict[wind[0]], 
                                            selected_dof_list, restraint_type)
            restrained_dof.plot_dof_vals(leg, run_nos, wind[0], wind[1], percent_traj_dict[wind[0]], 
                                            selected_dof_list, restraint_type)

        # RMSD for protein
        rmsd.plot_rmsds(leg, run_nos, percent_traj_dict, "protein")
        # RMSD for ligand
        rmsd.plot_rmsds(leg, run_nos, percent_traj_dict, "resname LIG and (not name H*)")
        # RMSD for syn-anti interconversion of ligand (ignore phenol group for which is rotatable and adds noise)
        rmsd.plot_rmsds(leg, run_nos, percent_traj_dict, "resname LIG and (not name H* OAA CAS CAD CAF CAG CAP CAE)")

    elif leg == "free":
        get_results.write_results(leg, run_nos)
        compare_conv.plot_stages_conv("analysis/convergence_data.pickle", leg)
        compare_conv.plot_overall_conv("analysis/convergence_data.pickle", leg)
        compare_pmfs.plot_all_pmfs(run_nos, leg)
        overlap.plot_overlap_mats(leg, run_nos)
        indiv_pmf_conv.plot_pmfs_conv(leg, run_nos)
        dh_dlam.plot_grads(leg, run_nos, percent_traj_dict, timestep, nrg_freq)

        # RMSD for syn-anti interconversion of ligand (ignore phenol group which is rotatable and adds noise)
        rmsd.plot_rmsds(leg, run_nos, percent_traj_dict, "resname LIG and (not name H* OAA CAS CAD CAF CAG CAP CAE)")


    print("###############################################################################################")
    print(f"Analysis of {leg} leg successfully completed!")