import numpy as np
import matplotlib.pyplot as plt

# Load data
data_2H = np.loadtxt('output/yield_phi_d_recon_data_2H.txt')[:,4]
sim_2H = np.loadtxt('output/yield_phi_d_recon_sim_2H.txt')[:,4]
tagged_2H = np.loadtxt('output/yield_phi_d_thrown_tagged_2H.txt')[:,4]

def lumi(energy_low, energy_high, target):
    if target == '2H':
        lumi_table = np.loadtxt('/work/halld2/home/boyu/src_analysis/flux/output/2H/lumi_summed_2H.txt')
    elif target == '4He':
        lumi_table = np.loadtxt('/work/halld2/home/boyu/src_analysis/flux/output/4He/lumi_summed_4He.txt')
    elif target == '12C':
        lumi_table = np.loadtxt('/work/halld2/home/boyu/src_analysis/flux/output/12C/lumi_summed_12C.txt')

    integrated_lumi = 0
    for i in range(len(lumi_table)):
        if (lumi_table[i,3] > energy_low) and (lumi_table[i,3] < energy_high):
            integrated_lumi += lumi_table[i][5]

    return integrated_lumi

minust_low = np.loadtxt('output/bin_edges.txt')[:,0]
minust_high = np.loadtxt('output/bin_edges.txt')[:,1]
energy_low = np.loadtxt('output/bin_edges.txt')[:,2]
energy_high = np.loadtxt('output/bin_edges.txt')[:,3]

acceptance_2H = np.zeros(len(minust_low))
dsdt_2H = np.zeros(len(minust_low))
err_2H = 1/np.sqrt(data_2H)

for i in range(len(minust_low)):
    if (sim_2H[i] == 0) or (tagged_2H[i] == 0):
        continue
    acceptance_2H[i] = sim_2H[i]/tagged_2H[i]
    dsdt_2H[i] = data_2H[i]/acceptance_2H[i]/lumi(energy_low[i], energy_high[i], '2H')/(minust_high[i]-minust_low[i])/1000

err_2H = dsdt_2H*err_2H

clas_t_low = np.array([0.350, 0.375, 0.400, 0.425, 0.450, 0.500, 0.550, 0.600, 0.700, 0.800, 1.000, 1.200, 1.400])
clas_t_high = np.array([0.375, 0.400, 0.425, 0.450, 0.500, 0.550, 0.600, 0.700, 0.800, 1.000, 1.200, 1.400, 2.000])
clas_cs_16 = np.array([10.21, 8.85, 7.32, 6.16, 4.73, 3.52, 2.66, 2.17, 1.40, 0.94, 0.57, 0.28, 0.19])
clas_stat_16 = np.array([0.82, 0.75, 0.59, 0.55, 0.34, 0.28, 0.24, 0.15, 0.12, 0.07, 0.06, 0.05, 0.02])
clas_syst_16 = np.array([1.70, 1.11, 0.94, 0.81, 0.60, 0.51, 0.38, 0.26, 0.16, 0.11, 0.07, 0.04, 0.03])
clas_cs_26 = np.array([8.63, 6.80, 4.57, 5.76, 3.99, 3.59, 2.11, 1.83, 1.32, 0.96, 0.57, 0.36, 0.15])
clas_stat_26 = np.array([0.80, 0.69, 0.53, 0.56, 0.33, 0.29, 0.22, 0.14, 0.12, 0.07, 0.05, 0.04, 0.02])
clas_syst_26 = np.array([1.04, 1.07, 0.74, 0.65, 0.55, 0.55, 0.28, 0.24, 0.20, 0.11, 0.06, 0.05, 0.02])

leps_t_low = np.array([0.38, 0.36, 0.34, 0.32, 0.30, 0.28, 0.26, 0.24, 0.22, 0.20, 0.18, 0.16, 0.14, 0.12, 0.10])
leps_t_high = np.array([0.4, 0.38, 0.36, 0.34, 0.32, 0.30, 0.28, 0.26, 0.24, 0.22, 0.20, 0.18, 0.16, 0.14, 0.12])
leps_cs = np.array([0.0005, 0.004, 0.0087, 0.0068, 0.0238, 0.0317, 0.0567, 0.0722, 0.092, 0.1186, 0.1749, 0.2033, 0.2544, 0.3101, 0.3396])*1000
leps_stat = np.array([0.0005, 0.002, 0.0035, 0.0028, 0.007, 0.0076, 0.0102, 0.0118, 0.0142, 0.0137, 0.0159, 0.0148, 0.0166, 0.0152, 0.0143])*1000

plt.errorbar((minust_low+minust_high)/2, data_2H, xerr=(minust_high-minust_low)/2, yerr=np.sqrt(data_2H), fmt='o')
plt.title(r"Yield of $d(\gamma, \phi d')$")
plt.xlabel(r'$-t[GeV^2/c]$')
plt.ylabel('Yield')
plt.savefig('output/fig_yield_2H.png')
plt.close()

plt.errorbar((minust_low+minust_high)/2, acceptance_2H, xerr=(minust_high-minust_low)/2, fmt='o')
plt.title(r"Acceptance of $d(\gamma, \phi d')$")
plt.xlabel(r'$-t[GeV^2/c]$')
plt.ylabel('Acceptance')
plt.savefig('output/fig_acceptance_2H.png')
plt.close()

plt.errorbar((minust_low+minust_high)/2, dsdt_2H, yerr=err_2H, fmt='.', label='SRC-CT, 8.2 GeV')
plt.legend()
plt.title(r"Differential cross section of $d(\gamma, \phi d')$")
plt.xlabel(r'$-t\:[GeV^2/c]$')
plt.ylabel(r'$d\sigma/dt\:[nb/(GeV^2/c)]$')
plt.yscale('log')
plt.savefig('output/fig_dsdt_2H_bare.png')
plt.close()

plt.errorbar((minust_low+minust_high)/2, dsdt_2H, yerr=err_2H, fmt='.', label='This work, SRC-CT, 8.2 GeV')
plt.errorbar((clas_t_low+clas_t_high)/2, clas_cs_16, yerr=clas_stat_16, fmt='.', label='CLAS 1.6-2.6 GeV')
plt.errorbar((clas_t_low+clas_t_high)/2, clas_cs_26, yerr=clas_stat_26, fmt='.', label='CLAS 2.6-3.6 GeV')
plt.errorbar((leps_t_low+leps_t_high)/2, leps_cs, yerr=leps_stat, fmt='.', label='LEPS 1.57-2.37 GeV')
plt.legend()
plt.title(r"Differential cross section of $d(\gamma, \phi d')$")
plt.xlabel(r'$-t\:[GeV^2/c]$')
plt.ylabel(r'$d\sigma/dt\:[nb/(GeV^2/c)]$')
plt.yscale('log')
plt.savefig('output/fig_dsdt_2H_compare.png')
plt.close()