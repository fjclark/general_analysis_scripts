import numpy as np
from math import pi, sin, log

# Constants
v0 = 1660.53907 # A^3, the standard state volume
kT = 0.593 # kcal mol-1 at 300 K

# Params
r = 5.92 #A
thetaA = 91.1 #deg
thetaB = 106 #deg

kr = 10 # kcal mol-1 A-2
kang = 10 # kcal mol-1 rad-2

# Calculation
numerator = 8*(pi**2)*v0*np.sqrt(kr*kang**5)
denominator = (r**2)*sin(thetaA*pi/180)*sin(thetaB*pi/180)*(2*pi*kT)**3

dg_r0 = -kT*log(numerator/denominator)

print(f"dg_r0 = {dg_r0:.2f} kcal mol-1")
print(numerator)
print(denominator)

