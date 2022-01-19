# Usage: check_error "directory containing .out files to check"

import sys
import os

def get_out_files(parent_dir):
    """Return list of .out files in parent dir"""
    out_files = []
    for file_name in os.listdir(parent_dir):
        if file_name.endswith(".out"):
            out_files.append(file_name)

    out_files.sort()

    return out_files

def get_lambdas(parent_dir):
    """Return list of lambda values based on files in parent dir"""
    lam_vals = []
    for file_name in os.listdir(parent_dir):
        if file_name.startswith("lambda"):
            lam_val = file_name.split(sep="-")[1]
            lam_vals.append(lam_val)
    
    lam_vals.sort()

    return lam_vals
    
def get_files_with_errors(file_list):
    """Return list of files from supplied list which terminate with error 255"""
    error_files = []
    for file in file_list:
        with open(file) as open_file:
            for line in open_file:
                if "exit code 255" in line:
                    error_files.append(file)
    return error_files


if __name__ == "__main__":
    parent_dir = sys.argv[1]
    out_files = get_out_files(parent_dir)
    lam_vals = get_lambdas(parent_dir)
    files_with_errors = get_files_with_errors(out_files)

    print("The lambda values at which errors occurred are:")

    for file in files_with_errors:
        index = int(file.split(sep=".")[1])
        print(lam_vals[index])