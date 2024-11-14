import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

mass_deuteron    = 1.875612

# Function to fit the forward cross section
def diff_cs_log(s, k, c):
    return k*s+c

def forward_cs(s, k, c, b):
    return k*np.power((s-b),c)

# t dependence
minus_t = (np.genfromtxt('input/theory_slope.txt', delimiter=','))[:,0]
log_sigma = (np.genfromtxt('input/theory_slope.txt', delimiter=','))[:,1]

plt.plot(minus_t, log_sigma, 'o')
popt_slope, pcov_slope = curve_fit(diff_cs_log, minus_t, log_sigma)
plt.plot(minus_t, diff_cs_log(minus_t, *popt_slope), 'r-')
np.savetxt("output/theory_slope.txt", popt_slope)
plt.savefig('output/theory_slope.png')
plt.close()

# Forward cross section
photon_energy = (np.genfromtxt('input/theory_forward.txt', delimiter=','))[:,0]
com_energy_sq = 2*photon_energy*mass_deuteron + mass_deuteron*mass_deuteron
sigma = (np.genfromtxt('input/theory_forward.txt', delimiter=','))[:,1]

plt.plot(com_energy_sq, sigma, 'o')
popt_forward, pcov_forward = curve_fit(forward_cs, com_energy_sq, sigma)
plt.plot(com_energy_sq, forward_cs(com_energy_sq, *popt_forward), 'r-')
np.savetxt("output/theory_forward.txt", popt_forward)
plt.savefig('output/theory_forward.png')
plt.close()

