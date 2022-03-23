"""Check that all calculations completed successfully; list failed windows if not"""

from .dir_paths import get_dir_paths
import glob

def check_success(run_nos = [1,2,3,4,5], leg="bound"):
    """Check that all lambda windows completed successfully; list
    failures if not.

    Args:
        run_nos (list, optional): Run numbers to check. Defaults to [1,2,3,4,5].
        leg (str, optional): bound or free. Defaults to bound.
    """
    failed_windows = []
    dir_paths = get_dir_paths(run_nos, leg)
    for run in dir_paths:
        for leg in dir_paths[run]:
            output_dir = dir_paths[run][leg]["output"]
            slurm_logs = glob.glob(f"{output_dir}/somd-array-gpu*")
            for log in slurm_logs:
                lam = "Undefined"
                with open(log, "rt") as f:
                    success = False
                    lines = f.readlines()
                    for l in lines:
                        if l.startswith("lambda is"):
                            lam = l.split()[-1]
                        if l.startswith("Simulation took"):
                            success = True
                            break
                    if not success:
                        failed_windows.append(f"{run} {leg} {lam}")

    if failed_windows == []:
        print("All windows ran successfully")
    else:
        print("The following windows failed:")
        for window in sorted(failed_windows):
            print(window)
        print(f"No failed windows = {len(failed_windows)}")
