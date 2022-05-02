# Calculates angle corresponding to energy penalty of 1kT at 298 K for
# force constant in kcal mol-1 rad-2
# Usage distance_at_kT "force constant"

import sys
from math import pi
from numpy import sqrt

# Constants
kT = 0.593 # kcal mol-1 at 298 K
k = float(sys.argv[1]) # kcal mol-1 rad-2
k_deg = k*(pi/180)**2 # convert force constant to kcal mol-1 deg-2

def angle_at_kT(force_const):
    """Calculate angle corresponding to energy penalty of 1kT at 298 K with
    supplied force constant"""
    return sqrt((2*kT)/force_const)
    
if __name__ == "__main__":
    print("Angle at 1kT at 298 K:")
    print(f"{angle_at_kT(k):.2f} rad")
    print(f"{angle_at_kT(k_deg):.0f} deg")
