import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

mass_deuteron    = 1.875612;

# Function to fit the differential cross section at a fixed photon energy
def func(t, a, b):
    return a * np.exp(b * t)

# Function to fit the forward cross section
def linear(s, k, c):
    return k*s+c

def geometric(s, k, c):
    return k*np.power(s,c)

def exponential(s, k, c):
    return k*np.exp(c*s)

# Function to fit the data at all energies
def func2(x, b, k, c):
    t, s = x
    return geometric(s, k, c) * np.exp(b * t)

# Read in the data points
energy = np.array([1.62, 1.72, 1.82, 1.92, 2.02, 2.12, 2.22, 2.32])
energy_com = 2*energy*mass_deuteron + mass_deuteron*mass_deuteron
energy_extended = np.arange(1.5, 10.5, 0.1)
energy_com_extended = 2*energy_extended*mass_deuteron + mass_deuteron*mass_deuteron
ds = []
ds.append(np.genfromtxt('input/Table1.csv', delimiter=','))
ds.append(np.genfromtxt('input/Table2.csv', delimiter=','))
ds.append(np.genfromtxt('input/Table3.csv', delimiter=','))
ds.append(np.genfromtxt('input/Table4.csv', delimiter=','))
ds.append(np.genfromtxt('input/Table5.csv', delimiter=','))
ds.append(np.genfromtxt('input/Table6.csv', delimiter=','))
ds.append(np.genfromtxt('input/Table7.csv', delimiter=','))
ds.append(np.genfromtxt('input/Table8.csv', delimiter=','))

# Check the ds/dt shape
plt.errorbar(ds[7][:,0], ds[7][:,3], yerr=ds[7][:,4], fmt='o')
plt.yscale('log')
plt.xlabel(r'$\tilde{t} (GeV^2)$')
plt.ylabel(r'$d\sigma/d\tilde{t} (\mu b/GeV^2)$')
plt.text(0.05, 0.9, r'$E_{\gamma}$='+str(energy[7]-0.05)+'-'+str(round(energy[7]+0.05, 2))+'GeV', fontsize=12, transform=plt.gca().transAxes)
plt.savefig('output/exponential_shape.png')
plt.close()

# Separate fit
popt = []
pcov = []
fig, ax = plt.subplots(2, 4, figsize=(20, 10), sharey=True, dpi=1000)
plt.subplots_adjust(hspace=0, wspace=0)
plt.yscale('log')
for i in range(2):
    ax[i,0].set_ylabel(r'$d\sigma/d\tilde{t} (\mu b/GeV^2)$')
for i in range(4):
    ax[1,i].set_xlabel(r'$\tilde{t} (GeV^2)$')
    ax[1,i].set_xticks(np.arange(-0.35,0.05,0.05))
    ax[0,i].set_xticks([])

for i in range(8):
    this_popt, this_pcov = curve_fit(func, ds[i][:,0], ds[i][:,3], sigma=ds[i][:,4], absolute_sigma=True)
    popt.append(this_popt)
    pcov.append(this_pcov)
    ax[i//4,i%4].errorbar(ds[i][:,0], ds[i][:,3], yerr=ds[i][:,4], fmt='o')
    ax[i//4,i%4].plot(ds[i][:,0], func(ds[i][:,0], *this_popt), 'r--')
    ax[i//4,i%4].set_xlim(-0.4,0)
    ax[i//4,i%4].text(0.05, 0.9, r'$E_{\gamma}$='+str(round(energy[i]-0.05, 2))+'-'+str(round(energy[i]+0.05, 2))+'GeV', fontsize=12, transform=ax[i//4,i%4].transAxes)
popt = np.array(popt)
plt.savefig('output/separate_fit.png')
plt.close()

plt.plot(energy_com, popt[:,1], 'o')
plt.plot(energy_com, np.average(popt[:,1])*np.ones(np.size(popt[:,1])), '-')
plt.xlabel(r'$s (GeV^2)$')
plt.ylabel(r'$b (GeV^{-2})$')
plt.ylim(0, 30)
plt.savefig('output/b_vs_E.png')
plt.close()

plt.plot(energy_com, popt[:,0], 'o')
popt_geometric, pcov_energy = curve_fit(geometric, energy_com, popt[:,0])
plt.plot(energy_com, geometric(energy_com, *popt_geometric), 'r-', label="geometric")
popt_linear, pcov_energy = curve_fit(linear, energy_com, popt[:,0])
plt.plot(energy_com, linear(energy_com, *popt_linear), 'g-', label="linear")
popt_exponential, pcov_energy = curve_fit(exponential, energy_com, popt[:,0])
plt.plot(energy_com, exponential(energy_com, *popt_exponential), 'b-', label="exponential")
plt.xlabel(r'$s (GeV^2)$')
plt.ylabel(r'$d\sigma/d\tilde{t} (\tilde{t}=0) (\mu b/GeV^2)$')
plt.xlim(energy_com[0], energy_com[-1])
plt.legend()
plt.savefig('output/dsdt0_vs_E.png')
plt.close()

plt.plot(energy_com, popt[:,0], 'o')
popt_geometric, pcov_energy = curve_fit(geometric, energy_com, popt[:,0])
plt.plot(energy_com_extended, geometric(energy_com_extended, *popt_geometric), 'r-', label="geometric")
popt_linear, pcov_energy = curve_fit(linear, energy_com, popt[:,0])
plt.plot(energy_com_extended, linear(energy_com_extended, *popt_linear), 'g-', label="linear")
popt_exponential, pcov_energy = curve_fit(exponential, energy_com, popt[:,0])
plt.plot(energy_com_extended, exponential(energy_com_extended, *popt_exponential), 'b-', label="exponential")
plt.xlabel(r'$s (GeV^2)$')
plt.ylabel(r'$d\sigma/d\tilde{t} (\tilde{t}=0) (\mu b/GeV^2)$')
plt.yscale('log')
plt.xlim(energy_com_extended[0], energy_com_extended[-1])
plt.legend()
plt.savefig('output/dsdt0_vs_E_extended.png')
plt.close()

# Combined fit
xdata = np.stack((ds[0][:,0], energy_com[0]*np.ones(ds[0][:,0].shape)), axis=1)
ydata = ds[0][:,3]
yerr = ds[0][:,4]
for i in range(1,8):
    xdata = np.concatenate((xdata, np.stack((ds[i][:,0], energy_com[i]*np.ones(ds[i][:,0].shape)), axis=1)), axis=0)
    ydata = np.concatenate((ydata, ds[i][:,3]))
    yerr = np.concatenate((yerr, ds[i][:,4]))
xdata = xdata.transpose()

popt_combined, pcov_combined = curve_fit(func2, xdata, ydata, sigma=yerr, absolute_sigma=True)

fig, ax = plt.subplots(2, 4, figsize=(20, 10), sharey=True, dpi=1000)
plt.subplots_adjust(hspace=0, wspace=0)
plt.yscale('log')
for i in range(2):
    ax[i,0].set_ylabel(r'$d\sigma/d\tilde{t} (\mu b/GeV^2)$')
for i in range(4):
    ax[1,i].set_xlabel(r'$\tilde{t} (GeV^2)$')
    ax[1,i].set_xticks(np.arange(-0.35,0.05,0.05))
    ax[0,i].set_xticks([])

for i in range(8):
    this_popt = np.array([geometric(energy_com[i], *popt_combined[1:]), popt_combined[0]])
    ax[i//4,i%4].errorbar(ds[i][:,0], ds[i][:,3], yerr=ds[i][:,4], fmt='o', label='data')
    ax[i//4,i%4].plot(ds[i][:,0], func(ds[i][:,0], *this_popt), 'r--', label='fit: a=%5.3f, b=%5.3f' % tuple(this_popt))
    ax[i//4,i%4].set_xlim(-0.4,0)
    ax[i//4,i%4].text(0.05, 0.9, r'$E_{\gamma}$='+str(round(energy_com[i]-0.05, 2))+'-'+str(round(energy_com[i]+0.05, 2))+'GeV', fontsize=12, transform=ax[i//4,i%4].transAxes)
plt.savefig('output/combined_fit.png')
plt.close()

print(popt_combined)