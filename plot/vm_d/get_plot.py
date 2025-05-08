import numpy as np
import matplotlib.pyplot as plt

def lumi(energy_min, energy_max, target):
    if target == '2H':
        lumi_table = np.loadtxt('/work/halld2/home/boyu/src_analysis/flux/output/2H/lumi_summed_2H.txt')
    elif target == '4He':
        lumi_table = np.loadtxt('/work/halld2/home/boyu/src_analysis/flux/output/4He/lumi_summed_4He.txt')
    elif target == '12C':
        lumi_table = np.loadtxt('/work/halld2/home/boyu/src_analysis/flux/output/12C/lumi_summed_12C.txt')

    integrated_lumi = 0
    for i in range(len(lumi_table)):
        if (lumi_table[i,3] > energy_min) and (lumi_table[i,3] < energy_max):
            integrated_lumi += lumi_table[i][5]

    return integrated_lumi

#======================================================================PHI_D_2H_DSDT======================================================================
# # Load data
# phi_d_2H_dsdt_yield_data = np.loadtxt('output/yield_phi_d_recon_data_2H_dsdt.txt')[:,4]
# phi_d_2H_dsdt_yield_sim = np.loadtxt('output/yield_phi_d_recon_sim_2H_dsdt.txt')[:,4]
# phi_d_2H_dsdt_yield_tagged = np.loadtxt('output/yield_phi_d_thrown_tagged_2H_dsdt.txt')[:,4]

# phi_d_2H_dsdt_energy_low = np.loadtxt('output/bins_phi_d_2H_dsdt.txt')[:,0]
# phi_d_2H_dsdt_energy_high = np.loadtxt('output/bins_phi_d_2H_dsdt.txt')[:,1]
# phi_d_2H_dsdt_minust_low = np.loadtxt('output/bins_phi_d_2H_dsdt.txt')[:,2]
# phi_d_2H_dsdt_minust_high = np.loadtxt('output/bins_phi_d_2H_dsdt.txt')[:,3]

# phi_d_2H_dsdt_acceptance = np.zeros(len(phi_d_2H_dsdt_minust_low))
# phi_d_2H_dsdt_result = np.zeros(len(phi_d_2H_dsdt_minust_low))
# phi_d_2H_dsdt_error_stat = 1/np.sqrt(phi_d_2H_dsdt_yield_data)

# for i in range(len(phi_d_2H_dsdt_minust_low)):
#     if (phi_d_2H_dsdt_yield_sim[i] == 0) or (phi_d_2H_dsdt_yield_tagged[i] == 0):
#         continue
#     phi_d_2H_dsdt_acceptance[i] = phi_d_2H_dsdt_yield_sim[i]/phi_d_2H_dsdt_yield_tagged[i]
#     phi_d_2H_dsdt_result[i] = phi_d_2H_dsdt_yield_data[i]/phi_d_2H_dsdt_acceptance[i]/lumi(phi_d_2H_dsdt_energy_low[i], phi_d_2H_dsdt_energy_high[i], '2H')/(phi_d_2H_dsdt_minust_high[i]-phi_d_2H_dsdt_minust_low[i])/0.489/1000

# phi_d_2H_dsdt_error_stat = phi_d_2H_dsdt_result*phi_d_2H_dsdt_error_stat

# clas_t_low = np.array([0.350, 0.375, 0.400, 0.425, 0.450, 0.500, 0.550, 0.600, 0.700, 0.800, 1.000, 1.200, 1.400])
# clas_t_high = np.array([0.375, 0.400, 0.425, 0.450, 0.500, 0.550, 0.600, 0.700, 0.800, 1.000, 1.200, 1.400, 2.000])
# clas_cs_16 = np.array([10.21, 8.85, 7.32, 6.16, 4.73, 3.52, 2.66, 2.17, 1.40, 0.94, 0.57, 0.28, 0.19])
# clas_stat_16 = np.array([0.82, 0.75, 0.59, 0.55, 0.34, 0.28, 0.24, 0.15, 0.12, 0.07, 0.06, 0.05, 0.02])
# clas_syst_16 = np.array([1.70, 1.11, 0.94, 0.81, 0.60, 0.51, 0.38, 0.26, 0.16, 0.11, 0.07, 0.04, 0.03])
# clas_cs_26 = np.array([8.63, 6.80, 4.57, 5.76, 3.99, 3.59, 2.11, 1.83, 1.32, 0.96, 0.57, 0.36, 0.15])
# clas_stat_26 = np.array([0.80, 0.69, 0.53, 0.56, 0.33, 0.29, 0.22, 0.14, 0.12, 0.07, 0.05, 0.04, 0.02])
# clas_syst_26 = np.array([1.04, 1.07, 0.74, 0.65, 0.55, 0.55, 0.28, 0.24, 0.20, 0.11, 0.06, 0.05, 0.02])

# leps_t_low = np.array([0.38, 0.36, 0.34, 0.32, 0.30, 0.28, 0.26, 0.24, 0.22, 0.20, 0.18, 0.16, 0.14, 0.12, 0.10])
# leps_t_high = np.array([0.4, 0.38, 0.36, 0.34, 0.32, 0.30, 0.28, 0.26, 0.24, 0.22, 0.20, 0.18, 0.16, 0.14, 0.12])
# leps_cs = np.array([0.0005, 0.004, 0.0087, 0.0068, 0.0238, 0.0317, 0.0567, 0.0722, 0.092, 0.1186, 0.1749, 0.2033, 0.2544, 0.3101, 0.3396])*1000
# leps_stat = np.array([0.0005, 0.002, 0.0035, 0.0028, 0.007, 0.0076, 0.0102, 0.0118, 0.0142, 0.0137, 0.0159, 0.0148, 0.0166, 0.0152, 0.0143])*1000

# plt.errorbar((phi_d_2H_dsdt_minust_low+phi_d_2H_dsdt_minust_high)/2, phi_d_2H_dsdt_yield_data, xerr=(phi_d_2H_dsdt_minust_high-phi_d_2H_dsdt_minust_low)/2, yerr=np.sqrt(phi_d_2H_dsdt_yield_data), fmt='o')
# plt.title(r"Yield of $d(\gamma, \phi d')$")
# plt.xlabel(r'$-t[GeV^2/c]$')
# plt.ylabel('Yield')
# plt.savefig('output/fig_phi_d_2H_dsdt_yield.png', dpi=300)
# plt.close()

# plt.errorbar((phi_d_2H_dsdt_minust_low+phi_d_2H_dsdt_minust_high)/2, phi_d_2H_dsdt_acceptance, xerr=(phi_d_2H_dsdt_minust_high-phi_d_2H_dsdt_minust_low)/2, fmt='o')
# plt.title(r"Acceptance of $d(\gamma, \phi d')$")
# plt.xlabel(r'$-t[GeV^2/c]$')
# plt.ylabel('Acceptance')
# plt.savefig('output/fig_phi_d_2H_dsdt_acceptance.png', dpi=300)
# plt.close()

# plt.errorbar((phi_d_2H_dsdt_minust_low+phi_d_2H_dsdt_minust_high)/2, phi_d_2H_dsdt_result, yerr=phi_d_2H_dsdt_error_stat, fmt='.', label='SRC-CT, 8.2 GeV')
# plt.legend()
# plt.title(r"Differential cross section of $d(\gamma, \phi d')$")
# plt.xlabel(r'$-t\ [GeV^2/c]$')
# plt.ylabel(r'$d\sigma/dt\ [nb/(GeV^2/c)]$')
# plt.yscale('log')
# plt.savefig('output/fig_phi_d_2H_dsdt_sigma_bare.png', dpi=300)
# plt.close()

# plt.errorbar((phi_d_2H_dsdt_minust_low+phi_d_2H_dsdt_minust_high)/2, phi_d_2H_dsdt_result, yerr=phi_d_2H_dsdt_error_stat, fmt='k.', label='This work, SRC-CT, 8.2 GeV')
# plt.errorbar((clas_t_low+clas_t_high)/2, clas_cs_16, yerr=clas_stat_16, fmt='s', markersize=4, fillstyle='none', label='CLAS 1.6-2.6 GeV')
# plt.errorbar((clas_t_low+clas_t_high)/2, clas_cs_26, yerr=clas_stat_26, fmt='s', markersize=4, fillstyle='none', label='CLAS 2.6-3.6 GeV')
# plt.errorbar((leps_t_low+leps_t_high)/2, leps_cs, yerr=leps_stat, fmt='s', markersize=4, fillstyle='none', label='LEPS 1.57-2.37 GeV')
# plt.legend()
# plt.title(r"Differential cross section of $d(\gamma, \phi d')$")
# plt.text(0.3, 0.15, 'preliminary', fontsize=15, color='r', style='italic', ha='center', va='center')
# plt.xlabel(r'$-t\ [GeV^2/c]$')
# plt.ylabel(r'$d\sigma/dt\ [nb/(GeV^2/c)]$')
# plt.yscale('log')
# plt.savefig('output/fig_phi_d_2H_dsdt_sigma_compare.png', dpi=300)
# plt.close()

#======================================================================PHI_D_2H_WCOSTHETA======================================================================
phi_d_2H_Wcostheta_yield_data = np.loadtxt('output/yield_phi_d_recon_data_2H_Wcostheta.txt')[:,6]
phi_d_2H_Wcostheta_yield_sim = np.loadtxt('output/yield_phi_d_recon_sim_2H_Wcostheta.txt')[:,6]
phi_d_2H_Wcostheta_yield_tagged = np.loadtxt('output/yield_phi_d_thrown_tagged_2H_Wcostheta.txt')[:,6]

phi_d_2H_Wcostheta_energy_low = np.loadtxt('output/bins_phi_d_2H_Wcostheta.txt')[:,0]
phi_d_2H_Wcostheta_energy_high = np.loadtxt('output/bins_phi_d_2H_Wcostheta.txt')[:,1]
phi_d_2H_Wcostheta_minust_low = np.loadtxt('output/bins_phi_d_2H_Wcostheta.txt')[:,2]
phi_d_2H_Wcostheta_minust_high = np.loadtxt('output/bins_phi_d_2H_Wcostheta.txt')[:,3]
phi_d_2H_Wcostheta_costheta_low = np.loadtxt('output/bins_phi_d_2H_Wcostheta.txt')[:,4]
phi_d_2H_Wcostheta_costheta_high = np.loadtxt('output/bins_phi_d_2H_Wcostheta.txt')[:,5]

phi_d_2H_Wcostheta_acceptance = np.zeros(len(phi_d_2H_Wcostheta_minust_low))
phi_d_2H_Wcostheta_result = np.zeros(len(phi_d_2H_Wcostheta_minust_low))
phi_d_2H_Wcostheta_error_stat = 1/np.sqrt(phi_d_2H_Wcostheta_yield_data)

for i in range(len(phi_d_2H_Wcostheta_minust_low)):
    if (phi_d_2H_Wcostheta_yield_sim[i] == 0) or (phi_d_2H_Wcostheta_yield_tagged[i] == 0):
        continue
    phi_d_2H_Wcostheta_acceptance[i] = phi_d_2H_Wcostheta_yield_sim[i]/phi_d_2H_Wcostheta_yield_tagged[i]
    phi_d_2H_Wcostheta_result[i] = phi_d_2H_Wcostheta_yield_data[i]/phi_d_2H_Wcostheta_acceptance[i]

phi_d_2H_Wcostheta_result[0:10] /= np.sum(phi_d_2H_Wcostheta_result[0:10])
phi_d_2H_Wcostheta_result[10:15] /= np.sum(phi_d_2H_Wcostheta_result[10:15])
phi_d_2H_Wcostheta_result /= phi_d_2H_Wcostheta_costheta_high - phi_d_2H_Wcostheta_costheta_low
phi_d_2H_Wcostheta_error_stat = phi_d_2H_Wcostheta_result*phi_d_2H_Wcostheta_error_stat

plt.figure(figsize=(15, 6))
plt.subplot(121)
plt.errorbar((phi_d_2H_Wcostheta_costheta_low[0:10]+phi_d_2H_Wcostheta_costheta_high[0:10])/2, phi_d_2H_Wcostheta_result[0:10], yerr=phi_d_2H_Wcostheta_error_stat[0:10], fmt='.', label='SRC-CT, 8.2 GeV')
plt.plot(np.linspace(-1, 1, 100), 0.75*(1-np.linspace(-1, 1, 100)**2), 'r--', label='SCHC+NPE')
plt.xlabel(r'$\cos\theta_{H}$')
plt.ylabel(r'$W(\cos\theta_{H})$')
plt.title(r'$0.2<-t<0.8\ \mathrm{GeV}^2$')
plt.subplot(122)
plt.errorbar((phi_d_2H_Wcostheta_costheta_low[10:15]+phi_d_2H_Wcostheta_costheta_high[10:15])/2, phi_d_2H_Wcostheta_result[10:15], yerr=phi_d_2H_Wcostheta_error_stat[10:15], fmt='.', label='SRC-CT, 8.2 GeV')
plt.plot(np.linspace(-1, 1, 100), 0.75*(1-np.linspace(-1, 1, 100)**2), 'r--', label='SCHC+NPE')
plt.xlabel(r'$\cos\theta_{H}$')
plt.ylabel(r'$W(\cos\theta_{H})$')
plt.title(r'$0.8<-t<2.0\ \mathrm{GeV}^2$')
plt.legend()
plt.savefig('output/fig_phi_d_2H_Wcostheta_result.png', dpi=300)
plt.close()

#======================================================================PHI_D_2H_Wphi======================================================================
phi_d_2H_Wphi_yield_data = np.loadtxt('output/yield_phi_d_recon_data_2H_Wphi.txt')[:,6]
phi_d_2H_Wphi_yield_sim = np.loadtxt('output/yield_phi_d_recon_sim_2H_Wphi.txt')[:,6]
phi_d_2H_Wphi_yield_tagged = np.loadtxt('output/yield_phi_d_thrown_tagged_2H_Wphi.txt')[:,6]

phi_d_2H_Wphi_energy_low = np.loadtxt('output/bins_phi_d_2H_Wphi.txt')[:,0]
phi_d_2H_Wphi_energy_high = np.loadtxt('output/bins_phi_d_2H_Wphi.txt')[:,1]
phi_d_2H_Wphi_minust_low = np.loadtxt('output/bins_phi_d_2H_Wphi.txt')[:,2]
phi_d_2H_Wphi_minust_high = np.loadtxt('output/bins_phi_d_2H_Wphi.txt')[:,3]
phi_d_2H_Wphi_phi_low = np.loadtxt('output/bins_phi_d_2H_Wphi.txt')[:,4]
phi_d_2H_Wphi_phi_high = np.loadtxt('output/bins_phi_d_2H_Wphi.txt')[:,5]

phi_d_2H_Wphi_acceptance = np.zeros(len(phi_d_2H_Wphi_minust_low))
phi_d_2H_Wphi_result = np.zeros(len(phi_d_2H_Wphi_minust_low))
phi_d_2H_Wphi_error_stat = 1/np.sqrt(phi_d_2H_Wphi_yield_data)

for i in range(len(phi_d_2H_Wphi_minust_low)):
    if (phi_d_2H_Wphi_yield_sim[i] == 0) or (phi_d_2H_Wphi_yield_tagged[i] == 0):
        continue
    phi_d_2H_Wphi_acceptance[i] = phi_d_2H_Wphi_yield_sim[i]/phi_d_2H_Wphi_yield_tagged[i]
    phi_d_2H_Wphi_result[i] = phi_d_2H_Wphi_yield_data[i]/phi_d_2H_Wphi_acceptance[i]

phi_d_2H_Wphi_result[0:9] /= np.sum(phi_d_2H_Wphi_result[0:9])
phi_d_2H_Wphi_result[9:18] /= np.sum(phi_d_2H_Wphi_result[9:18])
phi_d_2H_Wphi_result /= (phi_d_2H_Wphi_phi_high - phi_d_2H_Wphi_phi_low)/180*np.pi
phi_d_2H_Wphi_error_stat = phi_d_2H_Wphi_result*phi_d_2H_Wphi_error_stat

plt.figure(figsize=(15, 6))
plt.subplot(121)
plt.errorbar((phi_d_2H_Wphi_phi_low[0:9]+phi_d_2H_Wphi_phi_high[0:9])/2, 2*np.pi*phi_d_2H_Wphi_result[0:9], yerr=2*np.pi*phi_d_2H_Wphi_error_stat[0:9], fmt='.', label='SRC-CT, 8.2 GeV')
plt.plot(np.linspace(-180, 180, 360), np.ones(360), 'r--', label='SCHC+NPE')
plt.xlim(-180, 180)
plt.ylim(0, 2)
plt.xlabel(r'$\phi_{H}\ [\mathrm{deg}]$')
plt.ylabel(r'$2\pi W(\phi_{H})$')
plt.title(r'$0.2<-t<0.8\ \mathrm{GeV}^2$')
plt.subplot(122)
plt.errorbar((phi_d_2H_Wphi_phi_low[9:18]+phi_d_2H_Wphi_phi_high[9:18])/2, 2*np.pi*phi_d_2H_Wphi_result[9:18], yerr=2*np.pi*phi_d_2H_Wphi_error_stat[9:18], fmt='.', label='SRC-CT, 8.2 GeV')
plt.plot(np.linspace(-180, 180, 360), np.ones(360), 'r--', label='SCHC+NPE')
plt.xlim(-180, 180)
plt.ylim(0, 2)
plt.xlabel(r'$\phi_{H}\ [\mathrm{deg}]$')
plt.ylabel(r'$2\pi W(\phi_{H})$')
plt.title(r'$0.8<-t<2.0\ \mathrm{GeV}^2$')
plt.legend()
plt.savefig('output/fig_phi_d_2H_Wphi_result.png', dpi=300)
plt.close()