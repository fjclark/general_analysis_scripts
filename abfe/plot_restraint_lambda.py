import numpy as np
import matplotlib.pyplot as plt

k = 10 # kcal mol-1

def nrg(x,scale_fac=1):
    return 0.5*k*scale_fac*x**2

x = np.arange(-5,5,0.01)
scale_facts = np.arange(0,1.1,0.1)

fig, ax = plt.subplots()

ax.set_xlabel("Distance / $\AA$")
ax.set_ylabel("Energy / kcal mol$^{-1}$")

for scale_fac in scale_facts:
    y = list(map(lambda z: nrg(z, scale_fac), x))
    ax.plot(x,y,label = f"lambda: {scale_fac:.1f}")

ax.set_title(f"Scaling of Harmonic Distance Restraint with\nLambda for k = {k} kcal/ mol")
ax.set_ylim(0,2.5)
ax.legend()
ax.axhline(y=0.593, color='r', linestyle='--')
ax.text(-4.7, 0.65, "kT at 298 K", color='r')
plt.savefig("scaled_restraint_energies.png")
