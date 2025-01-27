import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

mass_deuteron    = 1.875612

# Function to fit the forward cross section
def diff_cs_log(s, k, c):
    return k*s+c

def forward_cs(s, k, c, b):
    return k*np.power(abs(s-b),c)

# t dependence
minus_t = (np.genfromtxt('input/theory_slope.txt', delimiter=','))[:,0]
log_sigma = (np.genfromtxt('input/theory_slope.txt', delimiter=','))[:,1]

plt.plot(minus_t, 10**(log_sigma), 'o', label="Titov calculation")
popt_slope, pcov_slope = curve_fit(diff_cs_log, minus_t, log_sigma)
plt.plot(minus_t, 10**(diff_cs_log(minus_t, *popt_slope)), 'r-', label="Parameterized fitting")
plt.xlabel(r"$-t (GeV^2)$")
plt.ylabel(r"$d\sigma /dt (\mu b/GeV^2)$")
plt.yscale('log')
plt.legend()
plt.title(r"$\gamma d \rightarrow \phi d, E_{\gamma}=2.07-2.17 \mathrm{ GeV}$")
np.savetxt("output/theory_slope.txt", popt_slope)
plt.savefig('output/theory_slope.png')
plt.close()

# Forward cross section
photon_energy = (np.genfromtxt('input/theory_forward.txt', delimiter=','))[:,0]
com_energy_sq = 2*photon_energy*mass_deuteron + mass_deuteron*mass_deuteron
sigma = (np.genfromtxt('input/theory_forward.txt', delimiter=','))[:,1]

fig = plt.figure(figsize=(10, 7))
ax1 = fig.add_subplot(111)
ax1.plot(com_energy_sq, sigma, 'o', label="Titov calculation")
popt_forward, pcov_forward = curve_fit(forward_cs, com_energy_sq, sigma)
ax1.plot(com_energy_sq, forward_cs(com_energy_sq, *popt_forward), 'r-', label="Parameterized fitting")
ax1.set_xlabel(r"$s (GeV^2)$")
ax1.set_ylabel(r"$d\sigma /dt (\mu b/GeV^2)$")
ax1.legend()
ax2 = ax1.twiny()
ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks(2*np.linspace(1.5, 2.5, 11)*mass_deuteron + mass_deuteron*mass_deuteron)
ax2.set_xticklabels(np.linspace(1.5, 2.5, 11))
ax2.set_xlabel(r"$E_{\gamma} (GeV)$")
plt.title(r"$\gamma d \rightarrow \phi d, \theta_{c.m.} = 0$")
np.savetxt("output/theory_forward.txt", popt_forward)
plt.savefig('output/theory_forward.png')
plt.close()