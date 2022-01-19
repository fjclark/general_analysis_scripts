# Calculates distance corresponding to energy penalty of 1kT at 298 K for
# force constant in kcal mol-1 A-2
# Usage distance_at_kT "force constant"

import sys
from numpy import sqrt

# Constants
kT = 0.593 # kcal mol-1 at 298 K
k = int(sys.argv[1])

def dist_at_kT(force_const):
    """Calculate distance (A) corresponding to energy penalty of 1kT at 298 K with
    supplied force constant"""
    return sqrt((2*kT)/force_const)

if __name__ == "__main__":
    print(f"Distance at 1kT and 298 K:\n{dist_at_kT(k):.2f} Angstrom")