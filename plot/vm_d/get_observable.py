import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages

rad_to_deg = 180/np.pi
file_pdf = PdfPages("/work/halld2/home/boyu/src_analysis/plot/vm_d/output/plots_vm_d.pdf")

def dsdt_func(minust, a1, b1, a2, b2):
    return a1*np.exp(-b1*minust) + a2*np.exp(-b2*minust)

def Wcostheta_func(costheta, c, alpha):
    return 0.75*((3*c-1)*costheta**2 + (1-c)) + alpha*costheta

def Wphi_func(phi, c):
    return 1-2*c*np.cos(2*phi*np.pi/180)

def lumi(energy_min, energy_max, target):
    if target == '2H':
        lumi_table = np.loadtxt('/work/halld2/home/boyu/src_analysis/flux/output/2H/lumi_summed_2H.txt')
    elif target == '4He':
        lumi_table = np.loadtxt('/work/halld2/home/boyu/src_analysis/flux/output/4He/lumi_summed_4He.txt')
    elif target == '12C':
        lumi_table = np.loadtxt('/work/halld2/home/boyu/src_analysis/flux/output/12C/lumi_summed_12C.txt')

    integrated_lumi = np.zeros(energy_min.shape, dtype=float)
    for i in range(len(energy_min)):
        for j in range(len(lumi_table)):
            if (lumi_table[j,3] > energy_min[i]) and (lumi_table[j,3] < energy_max[i]):
                integrated_lumi[i] += lumi_table[j][5]

    return integrated_lumi

def normalize_distribution(results, energy_bins, t_bins):
    for i in range(len(results)):
        if (i == 0):
            index = 0
        elif (i == len(results) - 1):
            results[index:i+1] /= np.sum(results[index:i+1])
        else:
            if (energy_bins[i] != energy_bins[i-1]) or (t_bins[i] != t_bins[i-1]):
                results[index:i] /= np.sum(results[index:i])
                index = i
    return results

#======================================================================phi_d_2H_dsdt======================================================================

# Read the bin edges
phi_d_2H_dsdt_energy_low            = np.loadtxt('configs/bins_phi_d_dsdt.txt')[:,0]
phi_d_2H_dsdt_energy_high           = np.loadtxt('configs/bins_phi_d_dsdt.txt')[:,1]
phi_d_2H_dsdt_energy_center         = np.loadtxt('output/yield_phi_d_recon_exc_data_2H_dsdt.txt')[:,0]
phi_d_2H_dsdt_energy_width          = np.loadtxt('output/yield_phi_d_recon_exc_data_2H_dsdt.txt')[:,1]
phi_d_2H_dsdt_energy_middle         = (phi_d_2H_dsdt_energy_high + phi_d_2H_dsdt_energy_low) / 2
phi_d_2H_dsdt_energy_size           = (phi_d_2H_dsdt_energy_high - phi_d_2H_dsdt_energy_low) / 2
phi_d_2H_dsdt_minust_low            = np.loadtxt('configs/bins_phi_d_dsdt.txt')[:,2]
phi_d_2H_dsdt_minust_high           = np.loadtxt('configs/bins_phi_d_dsdt.txt')[:,3]
phi_d_2H_dsdt_minust_center         = np.loadtxt('output/yield_phi_d_recon_exc_data_2H_dsdt.txt')[:,4]
phi_d_2H_dsdt_minust_width          = np.loadtxt('output/yield_phi_d_recon_exc_data_2H_dsdt.txt')[:,5]
phi_d_2H_dsdt_minust_middle         = (phi_d_2H_dsdt_minust_high + phi_d_2H_dsdt_minust_low) / 2
phi_d_2H_dsdt_minust_size           = (phi_d_2H_dsdt_minust_high - phi_d_2H_dsdt_minust_low) / 2

# Read the yield numbers
phi_d_2H_dsdt_yield_data            = np.loadtxt('output/yield_phi_d_recon_exc_data_2H_dsdt.txt')[:,8]
phi_d_2H_dsdt_yield_sim             = np.loadtxt('output/yield_phi_d_recon_exc_sim_2H_dsdt.txt')[:,8]
phi_d_2H_dsdt_yield_tagged          = np.loadtxt('output/yield_phi_d_thrown_exc_tagged_2H_dsdt.txt')[:,8]
phi_d_2H_dsdt_yield_data_statserr   = np.loadtxt('output/yield_phi_d_recon_exc_data_2H_dsdt.txt')[:,9]/phi_d_2H_dsdt_yield_data
phi_d_2H_dsdt_yield_sim_statserr    = np.loadtxt('output/yield_phi_d_recon_exc_sim_2H_dsdt.txt')[:,9]/phi_d_2H_dsdt_yield_sim
phi_d_2H_dsdt_yield_tagged_statserr = np.loadtxt('output/yield_phi_d_thrown_exc_tagged_2H_dsdt.txt')[:,9]/phi_d_2H_dsdt_yield_tagged

# Calculate the efficiency
phi_d_2H_dsdt_efficiency            = phi_d_2H_dsdt_yield_sim/phi_d_2H_dsdt_yield_tagged
phi_d_2H_dsdt_efficiency_statserr   = phi_d_2H_dsdt_efficiency*np.sqrt(phi_d_2H_dsdt_yield_sim_statserr**2 + phi_d_2H_dsdt_yield_tagged_statserr**2)

# Calculate the results
phi_d_2H_dsdt_results               = phi_d_2H_dsdt_yield_data/phi_d_2H_dsdt_efficiency/lumi(phi_d_2H_dsdt_energy_low, phi_d_2H_dsdt_energy_high, '2H')/(phi_d_2H_dsdt_minust_high-phi_d_2H_dsdt_minust_low)/0.489/1000
phi_d_2H_dsdt_results_statserr      = phi_d_2H_dsdt_results*np.sqrt(phi_d_2H_dsdt_yield_data_statserr**2 + phi_d_2H_dsdt_efficiency_statserr**2)

# Find the indices for the different energy and t bins
index = []
for i in range(len(phi_d_2H_dsdt_results)):
    if (i == 0):
        index.append(i)
    elif (i == len(phi_d_2H_dsdt_results) - 1):
        index.append(i+1)
    else:
        if (phi_d_2H_dsdt_energy_low[i] != phi_d_2H_dsdt_energy_low[i-1]):
            index.append(i)

print(index)

# Data points from CLAS
phi_d_2H_dsdt_clas_minust_low           = np.array([0.350, 0.375, 0.400, 0.425, 0.450, 0.500, 0.550, 0.600, 0.700, 0.800, 1.000, 1.200, 1.400])
phi_d_2H_dsdt_clas_minust_high          = np.array([0.375, 0.400, 0.425, 0.450, 0.500, 0.550, 0.600, 0.700, 0.800, 1.000, 1.200, 1.400, 2.000])
phi_d_2H_dsdt_clas_minust_center        = (phi_d_2H_dsdt_clas_minust_high + phi_d_2H_dsdt_clas_minust_low) / 2
phi_d_2H_dsdt_clas_minust_width         = (phi_d_2H_dsdt_clas_minust_high - phi_d_2H_dsdt_clas_minust_low) / 2
phi_d_2H_dsdt_clas_results_16           = np.array([10.21, 8.85, 7.32, 6.16, 4.73, 3.52, 2.66, 2.17, 1.40, 0.94, 0.57, 0.28, 0.19])
phi_d_2H_dsdt_clas_results_16_statserr  = np.array([0.82, 0.75, 0.59, 0.55, 0.34, 0.28, 0.24, 0.15, 0.12, 0.07, 0.06, 0.05, 0.02])
phi_d_2H_dsdt_clas_results_16_systerr   = np.array([1.70, 1.11, 0.94, 0.81, 0.60, 0.51, 0.38, 0.26, 0.16, 0.11, 0.07, 0.04, 0.03])
phi_d_2H_dsdt_clas_results_26           = np.array([8.63, 6.80, 4.57, 5.76, 3.99, 3.59, 2.11, 1.83, 1.32, 0.96, 0.57, 0.36, 0.15])
phi_d_2H_dsdt_clas_results_26_statserr  = np.array([0.80, 0.69, 0.53, 0.56, 0.33, 0.29, 0.22, 0.14, 0.12, 0.07, 0.05, 0.04, 0.02])
phi_d_2H_dsdt_clas_results_26_systerr   = np.array([1.04, 1.07, 0.74, 0.65, 0.55, 0.55, 0.28, 0.24, 0.20, 0.11, 0.06, 0.05, 0.02])

# Data points from LEPS
phi_d_2H_dsdt_leps_minust_low           = np.array([0.38, 0.36, 0.34, 0.32, 0.30, 0.28, 0.26, 0.24, 0.22, 0.20, 0.18, 0.16, 0.14, 0.12, 0.10])
phi_d_2H_dsdt_leps_minust_high          = np.array([0.4, 0.38, 0.36, 0.34, 0.32, 0.30, 0.28, 0.26, 0.24, 0.22, 0.20, 0.18, 0.16, 0.14, 0.12])
phi_d_2H_dsdt_leps_minust_center        = (phi_d_2H_dsdt_leps_minust_high + phi_d_2H_dsdt_leps_minust_low) / 2
phi_d_2H_dsdt_leps_minust_width         = (phi_d_2H_dsdt_leps_minust_high - phi_d_2H_dsdt_leps_minust_low) / 2
phi_d_2H_dsdt_leps_results_157          = np.array([0.0005, 0.004, 0.0087, 0.0068, 0.0238, 0.0317, 0.0567, 0.0722, 0.092, 0.1186, 0.1749, 0.2033, 0.2544, 0.3101, 0.3396])*1000
phi_d_2H_dsdt_leps_results_157_statserr = np.array([0.0005, 0.002, 0.0035, 0.0028, 0.007, 0.0076, 0.0102, 0.0118, 0.0142, 0.0137, 0.0159, 0.0148, 0.0166, 0.0152, 0.0143])*1000

# Theoretical prediction
# phi_d_2H_dsdt_theory_minust = np.loadtxt('ds_dt_theory/10mb_8gev_m0.csv', delimiter=',')[:,0]
# phi_d_2H_dsdt_theory_results_10mb_8gev_m0 = np.loadtxt('ds_dt_theory/10mb_8gev_m0.csv', delimiter=',')[:,1]
# phi_d_2H_dsdt_theory_results_10mb_8gev_m1 = np.loadtxt('ds_dt_theory/10mb_8gev_m1.csv', delimiter=',')[:,1]
# phi_d_2H_dsdt_theory_results_10mb_8gev = phi_d_2H_dsdt_theory_results_10mb_8gev_m1*2/3 + phi_d_2H_dsdt_theory_results_10mb_8gev_m0/3
# phi_d_2H_dsdt_theory_results_30mb_8gev_m0 = np.loadtxt('ds_dt_theory/30mb_8gev_m0.csv', delimiter=',')[:,1]
# phi_d_2H_dsdt_theory_results_30mb_8gev_m1 = np.loadtxt('ds_dt_theory/30mb_8gev_m1.csv', delimiter=',')[:,1]
# phi_d_2H_dsdt_theory_results_30mb_8gev = phi_d_2H_dsdt_theory_results_30mb_8gev_m1*2/3 + phi_d_2H_dsdt_theory_results_30mb_8gev_m0/3

# Plot the data yield
fig = plt.figure(figsize=(8, 6), dpi=300)
plt.errorbar(phi_d_2H_dsdt_minust_center[index[0]:index[1]], phi_d_2H_dsdt_yield_data[index[0]:index[1]], xerr=phi_d_2H_dsdt_minust_width[index[0]:index[1]], yerr=phi_d_2H_dsdt_yield_data_statserr[index[0]:index[1]], fmt='b.', label='6-8 GeV')
plt.errorbar(phi_d_2H_dsdt_minust_center[index[1]:index[2]], phi_d_2H_dsdt_yield_data[index[1]:index[2]], xerr=phi_d_2H_dsdt_minust_width[index[1]:index[2]], yerr=phi_d_2H_dsdt_yield_data_statserr[index[1]:index[2]], fmt='k.', label='8-9 GeV')
plt.errorbar(phi_d_2H_dsdt_minust_center[index[2]:index[3]], phi_d_2H_dsdt_yield_data[index[2]:index[3]], xerr=phi_d_2H_dsdt_minust_width[index[2]:index[3]], yerr=phi_d_2H_dsdt_yield_data_statserr[index[2]:index[3]], fmt='r.', label='9-11 GeV')
plt.title(r"$d(\gamma, \phi d')$ yield vs $-t$")
plt.xlabel(r'$-t[GeV^2/c]$')
plt.ylabel(r'$\mathrm{Yield}$')
plt.xlim(0, 2)
plt.ylim(0, 200)
plt.legend()
file_pdf.savefig()
plt.close()

# Plot the efficiency
fig = plt.figure(figsize=(8, 6), dpi=300)
plt.errorbar(phi_d_2H_dsdt_minust_center, phi_d_2H_dsdt_efficiency, xerr=phi_d_2H_dsdt_minust_width, yerr=phi_d_2H_dsdt_efficiency_statserr, fmt='k.', label='Weighted')
plt.title(r"$d(\gamma, \phi d')$ efficiency vs $-t$")
plt.xlabel(r'$-t[GeV^2/c]$')
plt.ylabel(r'$\mathrm{Efficiency}$')
plt.xlim(0, 2)
plt.ylim(0, 0.4)
plt.legend()
file_pdf.savefig()
plt.close()

# Plot the results
fig = plt.figure(figsize=(8, 6), dpi=300)
plt.errorbar(phi_d_2H_dsdt_minust_center[index[0]:index[1]],       phi_d_2H_dsdt_results[index[0]:index[1]],          xerr=phi_d_2H_dsdt_minust_width[index[0]:index[1]],        yerr=phi_d_2H_dsdt_results_statserr[index[0]:index[1]],            fmt='b.', label='6-8 GeV')
plt.errorbar(phi_d_2H_dsdt_minust_center[index[1]:index[2]],       phi_d_2H_dsdt_results[index[1]:index[2]],          xerr=phi_d_2H_dsdt_minust_width[index[1]:index[2]],        yerr=phi_d_2H_dsdt_results_statserr[index[1]:index[2]],            fmt='k.', label='8-9 GeV')
plt.errorbar(phi_d_2H_dsdt_minust_center[index[2]:index[3]],       phi_d_2H_dsdt_results[index[2]:index[3]],          xerr=phi_d_2H_dsdt_minust_width[index[2]:index[3]],        yerr=phi_d_2H_dsdt_results_statserr[index[2]:index[3]],            fmt='r.', label='9-11 GeV')
plt.errorbar(phi_d_2H_dsdt_clas_minust_center,  phi_d_2H_dsdt_clas_results_16,  xerr=phi_d_2H_dsdt_clas_minust_width,   yerr=phi_d_2H_dsdt_clas_results_16_statserr,    fmt='ys', markersize=4, fillstyle='none', label='CLAS 1.6-2.6 GeV')
plt.errorbar(phi_d_2H_dsdt_clas_minust_center,  phi_d_2H_dsdt_clas_results_26,  xerr=phi_d_2H_dsdt_clas_minust_width,   yerr=phi_d_2H_dsdt_clas_results_26_statserr,    fmt='cs', markersize=4, fillstyle='none', label='CLAS 2.6-3.6 GeV')
plt.errorbar(phi_d_2H_dsdt_leps_minust_center,  phi_d_2H_dsdt_leps_results_157, xerr=phi_d_2H_dsdt_leps_minust_width,   yerr=phi_d_2H_dsdt_leps_results_157_statserr,   fmt='gs', markersize=4, fillstyle='none', label='LEPS 1.57-2.37 GeV')
# curve_fit_params, curve_fit_cov = curve_fit(dsdt_func, phi_d_2H_dsdt_minust_center, phi_d_2H_dsdt_results, p0=[0, 0, 0, 0])
# curve_fit_residuals = phi_d_2H_dsdt_results - dsdt_func(phi_d_2H_dsdt_minust_center, curve_fit_params[0], curve_fit_params[1], curve_fit_params[2], curve_fit_params[3])
# reduced_chi2 = np.sum((curve_fit_residuals/phi_d_2H_dsdt_results_statserr)**2)/(len(phi_d_2H_dsdt_results)-4)
# plt.plot(np.linspace(0, 2, 100), dsdt_func(np.linspace(0, 2, 100), curve_fit_params[0], curve_fit_params[1], curve_fit_params[2], curve_fit_params[3]), 'b--', label='Fit')
# plt.plot(phi_d_2H_dsdt_theory_minust, phi_d_2H_dsdt_theory_results_10mb_8gev, 'r--', label='10mb prediction')
# plt.plot(phi_d_2H_dsdt_theory_minust, phi_d_2H_dsdt_theory_results_30mb_8gev, 'b--', label='30mb prediction')
plt.text(0.3, 0.15, 'preliminary', fontsize=15, color='r', style='italic', ha='center', va='center')
plt.title(r"$d(\gamma, \phi d')$ differential cross section vs $-t$")
plt.xlabel(r'$-t\ [GeV^2/c]$')
plt.ylabel(r'$d\sigma/dt\ [nb/(GeV^2/c)]$')
plt.xlim(0, 2)
plt.ylim(1e-1, 1e3)
plt.yscale('log')
plt.legend()
file_pdf.savefig()
plt.close()

fig = plt.figure(figsize=(8, 6), dpi=300)
color_code = ['b', 'k', 'r']
plt.errorbar(phi_d_2H_dsdt_minust_center[index[0]:index[1]],       phi_d_2H_dsdt_results[index[0]:index[1]],          xerr=phi_d_2H_dsdt_minust_width[index[0]:index[1]],        yerr=phi_d_2H_dsdt_results_statserr[index[0]:index[1]],            fmt='b.', label='6-8 GeV')
plt.errorbar(phi_d_2H_dsdt_minust_center[index[1]:index[2]],       phi_d_2H_dsdt_results[index[1]:index[2]],          xerr=phi_d_2H_dsdt_minust_width[index[1]:index[2]],        yerr=phi_d_2H_dsdt_results_statserr[index[1]:index[2]],            fmt='k.', label='8-9 GeV')
plt.errorbar(phi_d_2H_dsdt_minust_center[index[2]:index[3]],       phi_d_2H_dsdt_results[index[2]:index[3]],          xerr=phi_d_2H_dsdt_minust_width[index[2]:index[3]],        yerr=phi_d_2H_dsdt_results_statserr[index[2]:index[3]],            fmt='r.', label='9-11 GeV')
plt.text(0.3, 0.15, 'preliminary', fontsize=15, color='r', style='italic', ha='center', va='center')
plt.title(r"$d(\gamma, \phi d')$ differential cross section vs $-t$")
plt.xlabel(r'$-t\ [GeV^2/c]$')
plt.ylabel(r'$d\sigma/dt\ [nb/(GeV^2/c)]$')
plt.xlim(0, 2)
plt.ylim(1e-1, 1e3)
plt.yscale('log')
plt.legend()
file_pdf.savefig()
plt.close()

fig = plt.figure(figsize=(8, 6), dpi=300)
color_code = ['b', 'k', 'r']
plt.errorbar(phi_d_2H_dsdt_minust_center[index[0]:index[1]],       phi_d_2H_dsdt_results[index[0]:index[1]],          xerr=phi_d_2H_dsdt_minust_width[index[0]:index[1]],        yerr=phi_d_2H_dsdt_results_statserr[index[0]:index[1]],            fmt='b.', label='6-8 GeV')
plt.errorbar(phi_d_2H_dsdt_minust_center[index[1]:index[2]],       phi_d_2H_dsdt_results[index[1]:index[2]],          xerr=phi_d_2H_dsdt_minust_width[index[1]:index[2]],        yerr=phi_d_2H_dsdt_results_statserr[index[1]:index[2]],            fmt='k.', label='8-9 GeV')
plt.errorbar(phi_d_2H_dsdt_minust_center[index[2]:index[3]],       phi_d_2H_dsdt_results[index[2]:index[3]],          xerr=phi_d_2H_dsdt_minust_width[index[2]:index[3]],        yerr=phi_d_2H_dsdt_results_statserr[index[2]:index[3]],            fmt='r.', label='9-11 GeV')
for i in range(3):
    curve_fit_params, curve_fit_cov = curve_fit(dsdt_func, phi_d_2H_dsdt_minust_center[index[i]:index[i+1]], phi_d_2H_dsdt_results[index[i]:index[i+1]], p0=[0, 0, 0, 0])
    curve_fit_residuals = phi_d_2H_dsdt_results[index[i]:index[i+1]] - dsdt_func(phi_d_2H_dsdt_minust_center[index[i]:index[i+1]], curve_fit_params[0], curve_fit_params[1], curve_fit_params[2], curve_fit_params[3])
    reduced_chi2 = np.sum((curve_fit_residuals/phi_d_2H_dsdt_results_statserr[index[i]:index[i+1]])**2)/(len(phi_d_2H_dsdt_results[index[i]:index[i+1]])-4)
    plt.plot(np.linspace(0, 2, 100), dsdt_func(np.linspace(0, 2, 100), curve_fit_params[0], curve_fit_params[1], curve_fit_params[2], curve_fit_params[3]), '--', color = color_code[i], label='Paras: %.2f, %.2f, %.2f, %.2f, $\chi^2$/ndf = %.2f' % (curve_fit_params[0], curve_fit_params[1], curve_fit_params[2], curve_fit_params[3], reduced_chi2))
plt.text(0.3, 0.15, 'preliminary', fontsize=15, color='r', style='italic', ha='center', va='center')
plt.title(r"$d(\gamma, \phi d')$ differential cross section vs $-t$")
plt.xlabel(r'$-t\ [GeV^2/c]$')
plt.ylabel(r'$d\sigma/dt\ [nb/(GeV^2/c)]$')
plt.xlim(0, 2)
plt.ylim(1e-1, 1e3)
plt.yscale('log')
plt.legend()
file_pdf.savefig()
plt.close()

fig = plt.figure(figsize=(8, 6), dpi=300)
color_code = ['b', 'k', 'r']
plt.errorbar(phi_d_2H_dsdt_minust_center[index[0]:index[1]],       phi_d_2H_dsdt_results[index[0]:index[1]],          xerr=phi_d_2H_dsdt_minust_width[index[0]:index[1]],        yerr=phi_d_2H_dsdt_results_statserr[index[0]:index[1]],            fmt='b.', label='6-8 GeV')
plt.errorbar(phi_d_2H_dsdt_minust_center[index[1]:index[2]],       10*phi_d_2H_dsdt_results[index[1]:index[2]],          xerr=phi_d_2H_dsdt_minust_width[index[1]:index[2]],        yerr=10*phi_d_2H_dsdt_results_statserr[index[1]:index[2]],            fmt='k.', label='8-9 GeV')
plt.errorbar(phi_d_2H_dsdt_minust_center[index[2]:index[3]],       100*phi_d_2H_dsdt_results[index[2]:index[3]],          xerr=phi_d_2H_dsdt_minust_width[index[2]:index[3]],        yerr=100*phi_d_2H_dsdt_results_statserr[index[2]:index[3]],            fmt='r.', label='9-11 GeV')
for i in range(3):
    curve_fit_params, curve_fit_cov = curve_fit(dsdt_func, phi_d_2H_dsdt_minust_center[index[i]:index[i+1]], phi_d_2H_dsdt_results[index[i]:index[i+1]], p0=[0, 0, 0, 0])
    curve_fit_residuals = phi_d_2H_dsdt_results[index[i]:index[i+1]] - dsdt_func(phi_d_2H_dsdt_minust_center[index[i]:index[i+1]], curve_fit_params[0], curve_fit_params[1], curve_fit_params[2], curve_fit_params[3])
    reduced_chi2 = np.sum((curve_fit_residuals/phi_d_2H_dsdt_results_statserr[index[i]:index[i+1]])**2)/(len(phi_d_2H_dsdt_results[index[i]:index[i+1]])-4)
    plt.plot(np.linspace(0, 2, 100), pow(10,i)*dsdt_func(np.linspace(0, 2, 100), curve_fit_params[0], curve_fit_params[1], curve_fit_params[2], curve_fit_params[3]), '--', color = color_code[i])
plt.text(0.3, 0.15, 'preliminary', fontsize=15, color='r', style='italic', ha='center', va='center')
plt.title(r"$d(\gamma, \phi d')$ differential cross section vs $-t$")
plt.xlabel(r'$-t\ [GeV^2/c]$')
plt.ylabel(r'$d\sigma/dt\ [nb/(GeV^2/c)]$')
plt.xlim(0, 2)
plt.ylim(1e-1, 1e4)
plt.yscale('log')
plt.legend()
file_pdf.savefig()
plt.close()

# #======================================================================phi_d_2H_Wcostheta======================================================================

# # Read the bin edges
# phi_d_2H_Wcostheta_energy_low               = np.loadtxt('configs/bins_phi_d_2H_Wcostheta.txt')[:,0]
# phi_d_2H_Wcostheta_energy_high              = np.loadtxt('configs/bins_phi_d_2H_Wcostheta.txt')[:,1]
# phi_d_2H_Wcostheta_minust_low               = np.loadtxt('configs/bins_phi_d_2H_Wcostheta.txt')[:,2]
# phi_d_2H_Wcostheta_minust_high              = np.loadtxt('configs/bins_phi_d_2H_Wcostheta.txt')[:,3]
# phi_d_2H_Wcostheta_costheta_low             = np.loadtxt('configs/bins_phi_d_2H_Wcostheta.txt')[:,4]
# phi_d_2H_Wcostheta_costheta_high            = np.loadtxt('configs/bins_phi_d_2H_Wcostheta.txt')[:,5]
# phi_d_2H_Wcostheta_costheta_center          = (phi_d_2H_Wcostheta_costheta_high + phi_d_2H_Wcostheta_costheta_low) / 2
# phi_d_2H_Wcostheta_costheta_width           = (phi_d_2H_Wcostheta_costheta_high - phi_d_2H_Wcostheta_costheta_low) / 2

# # Read the yield numbers
# phi_d_2H_Wcostheta_yield_data               = np.loadtxt('output/yield_phi_d_recon_data_2H_exc_Wcostheta.txt')[:,6]
# phi_d_2H_Wcostheta_yield_sim                = np.loadtxt('output/yield_phi_d_recon_sim_2H_exc_Wcostheta.txt')[:,6]
# phi_d_2H_Wcostheta_yield_tagged             = np.loadtxt('output/yield_phi_d_thrown_tagged_2H_Wcostheta.txt')[:,6]
# phi_d_2H_Wcostheta_yield_data_statserr      = np.sqrt(phi_d_2H_Wcostheta_yield_data)
# phi_d_2H_Wcostheta_yield_sim_statserr       = np.sqrt(phi_d_2H_Wcostheta_yield_sim)
# phi_d_2H_Wcostheta_yield_tagged_statserr    = np.sqrt(phi_d_2H_Wcostheta_yield_tagged)

# # Calculate the efficiency
# phi_d_2H_Wcostheta_efficiency               = phi_d_2H_Wcostheta_yield_sim/phi_d_2H_Wcostheta_yield_tagged
# phi_d_2H_Wcostheta_efficiency_statserr      = phi_d_2H_Wcostheta_efficiency * np.sqrt(1/phi_d_2H_Wcostheta_yield_sim + 1/phi_d_2H_Wcostheta_yield_tagged)

# # Calculate the results
# phi_d_2H_Wcostheta_results                  = phi_d_2H_Wcostheta_yield_data/phi_d_2H_Wcostheta_efficiency  # raw results
# phi_d_2H_Wcostheta_results                  = normalize_distribution(phi_d_2H_Wcostheta_results, phi_d_2H_Wcostheta_energy_low, phi_d_2H_Wcostheta_minust_low) # normalize to have the sum equal to 1
# phi_d_2H_Wcostheta_results                  = phi_d_2H_Wcostheta_results/(phi_d_2H_Wcostheta_costheta_high - phi_d_2H_Wcostheta_costheta_low)  # normalize to have the integral equal to 1
# phi_d_2H_Wcostheta_results_statserr         = phi_d_2H_Wcostheta_results/np.sqrt(phi_d_2H_Wcostheta_yield_data)

# # Find the indices for the different energy and t bins
# index = []
# for i in range(len(phi_d_2H_Wcostheta_results)):
#     if (i == 0):
#         index.append(i)
#     elif (i == len(phi_d_2H_Wcostheta_results) - 1):
#         index.append(i+1)
#     else:
#         if (phi_d_2H_Wcostheta_energy_low[i] != phi_d_2H_Wcostheta_energy_low[i-1]) or (phi_d_2H_Wcostheta_minust_low[i] != phi_d_2H_Wcostheta_minust_low[i-1]):
#             index.append(i)

# # Plot the data yield
# fig = plt.figure(figsize=(6*(len(index) - 1), 6))
# gs = fig.add_gridspec(1, len(index) - 1, wspace=0)
# axs = gs.subplots(sharex=True, sharey=True)
# for i in range(len(index) - 1):
#     axs[i].errorbar(phi_d_2H_Wcostheta_costheta_center[index[i]:index[i+1]], phi_d_2H_Wcostheta_yield_data[index[i]:index[i+1]], xerr=phi_d_2H_Wcostheta_costheta_width[index[i]:index[i+1]], yerr=phi_d_2H_Wcostheta_yield_data_statserr[index[i]:index[i+1]], fmt='k.', label='This work')
#     axs[i].set_title(r'$%.1f<\mathrm{E}_{\gamma}<%.1f\ \mathrm{GeV},\  %.1f<-t<%.1f\ \mathrm{GeV}^2$' % (phi_d_2H_Wcostheta_energy_low[index[i]], phi_d_2H_Wcostheta_energy_high[index[i]], phi_d_2H_Wcostheta_minust_low[index[i]], phi_d_2H_Wcostheta_minust_high[index[i]]))
# axs[0].set_xlim(-1, 1)
# axs[0].set_ylim(0, 150)
# axs[0].set_xticks(np.arange(-0.75, 0.9, 0.25))
# axs[0].set_ylabel(r'$\mathrm{Yield}$')
# fig.suptitle(r"$d(\gamma, \phi d')$ yield vs $\cos\vartheta$")
# fig.supxlabel(r'$\cos\vartheta$')
# plt.legend()
# plt.savefig('output/fig_phi_d_2H_Wcostheta_yield.png', dpi=300)
# plt.close()

# # Plot the efficiency
# fig = plt.figure(figsize=(6*(len(index) - 1), 6))
# gs = fig.add_gridspec(1, len(index) - 1, wspace=0)
# axs = gs.subplots(sharex=True, sharey=True)
# for i in range(len(index) - 1):
#     axs[i].errorbar(phi_d_2H_Wcostheta_costheta_center[index[i]:index[i+1]], phi_d_2H_Wcostheta_efficiency[index[i]:index[i+1]], xerr=phi_d_2H_Wcostheta_costheta_width[index[i]:index[i+1]], yerr=phi_d_2H_Wcostheta_efficiency_statserr[index[i]:index[i+1]], fmt='k.', label='This work')
#     axs[i].set_title(r'$%.1f<\mathrm{E}_{\gamma}<%.1f\ \mathrm{GeV},\  %.1f<-t<%.1f\ \mathrm{GeV}^2$' % (phi_d_2H_Wcostheta_energy_low[index[i]], phi_d_2H_Wcostheta_energy_high[index[i]], phi_d_2H_Wcostheta_minust_low[index[i]], phi_d_2H_Wcostheta_minust_high[index[i]]))
# axs[0].set_xlim(-1, 1)
# axs[0].set_ylim(0, 0.1)
# axs[0].set_xticks(np.arange(-0.75, 0.9, 0.25))
# axs[0].set_ylabel(r'$\mathrm{Efficiency}$')
# fig.suptitle(r"$d(\gamma, \phi d')$ efficiency vs $\cos\vartheta$")
# fig.supxlabel(r'$\cos\vartheta$')
# plt.legend()
# plt.savefig('output/fig_phi_d_2H_Wcostheta_efficiency.png', dpi=300)
# plt.close()

# # Plot the results
# fig = plt.figure(figsize=(6*(len(index) - 1), 6))
# gs = fig.add_gridspec(1, len(index) - 1, wspace=0)
# axs = gs.subplots(sharex=True, sharey=True)
# for i in range(len(index) - 1):
#     axs[i].errorbar(phi_d_2H_Wcostheta_costheta_center[index[i]:index[i+1]], phi_d_2H_Wcostheta_results[index[i]:index[i+1]], xerr=phi_d_2H_Wcostheta_costheta_width[index[i]:index[i+1]], yerr=phi_d_2H_Wcostheta_results_statserr[index[i]:index[i+1]], fmt='k.', label='This work')
#     axs[i].set_title(r'$%.1f<\mathrm{E}_{\gamma}<%.1f\ \mathrm{GeV},\  %.1f<-t<%.1f\ \mathrm{GeV}^2$' % (phi_d_2H_Wcostheta_energy_low[index[i]], phi_d_2H_Wcostheta_energy_high[index[i]], phi_d_2H_Wcostheta_minust_low[index[i]], phi_d_2H_Wcostheta_minust_high[index[i]]))
#     curve_fit_params, curve_fit_cov = curve_fit(Wcostheta_func, phi_d_2H_Wcostheta_costheta_center[index[i]:index[i+1]], phi_d_2H_Wcostheta_results[index[i]:index[i+1]], p0=[0.0, 0.0])
#     curve_fit_residuals = phi_d_2H_Wcostheta_results[index[i]:index[i+1]] - Wcostheta_func(phi_d_2H_Wcostheta_costheta_center[index[i]:index[i+1]], curve_fit_params[0], curve_fit_params[1])
#     reduced_chi2 = np.sum((curve_fit_residuals/phi_d_2H_Wcostheta_results_statserr[index[i]:index[i+1]])**2)/(len(phi_d_2H_Wcostheta_results[index[i]:index[i+1]])-2)
#     axs[i].plot(np.linspace(-1, 1, 100), Wcostheta_func(np.linspace(-1, 1, 100), curve_fit_params[0], curve_fit_params[1]), 'b--', label='Fit')
#     axs[i].text(-0.9, 0.95, r'$\rho^0_{00}=%.2f\pm%.2f$' % (curve_fit_params[0], np.sqrt(curve_fit_cov[0][0])), fontsize=10, color='b', ha='left', va='top')
#     axs[i].text(-0.9, 0.90, r'$\alpha=%.2e\pm%.2e$' % (curve_fit_params[1], np.sqrt(curve_fit_cov[1][1])), fontsize=10, color='b', ha='left', va='top')
#     axs[i].text(-0.9, 0.85, r'$\chi^2/\mathrm{dof}=%.2f$' % reduced_chi2, fontsize=10, color='b', ha='left', va='top')
#     axs[i].plot(np.linspace(-1, 1, 100), 0.75*(1-np.linspace(-1, 1, 100)**2), 'r--', label='SCHC+NPE')
# axs[0].set_xlim(-1, 1)
# axs[0].set_ylim(0, 1)
# axs[0].set_xticks(np.arange(-0.75, 0.9, 0.25))
# axs[0].set_ylabel(r'$W(\cos\vartheta)$')
# fig.suptitle(r"$d(\gamma, \phi d')$ normalized distribution of $\cos\vartheta$$")
# fig.supxlabel(r'$\cos\vartheta$')
# plt.legend()
# plt.savefig('output/fig_phi_d_2H_Wcostheta_results.png', dpi=300)
# plt.close()

# #======================================================================phi_d_2H_Wphi======================================================================

# # Read the bin edges
# phi_d_2H_Wphi_energy_low            = np.loadtxt('configs/bins_phi_d_2H_Wphi.txt')[:,0]
# phi_d_2H_Wphi_energy_high           = np.loadtxt('configs/bins_phi_d_2H_Wphi.txt')[:,1]
# phi_d_2H_Wphi_minust_low            = np.loadtxt('configs/bins_phi_d_2H_Wphi.txt')[:,2]
# phi_d_2H_Wphi_minust_high           = np.loadtxt('configs/bins_phi_d_2H_Wphi.txt')[:,3]
# phi_d_2H_Wphi_phi_low               = np.loadtxt('configs/bins_phi_d_2H_Wphi.txt')[:,4]
# phi_d_2H_Wphi_phi_high              = np.loadtxt('configs/bins_phi_d_2H_Wphi.txt')[:,5]
# phi_d_2H_Wphi_phi_center            = (phi_d_2H_Wphi_phi_high + phi_d_2H_Wphi_phi_low) / 2
# phi_d_2H_Wphi_phi_width             = (phi_d_2H_Wphi_phi_high - phi_d_2H_Wphi_phi_low) / 2

# # Read the yield numbers
# phi_d_2H_Wphi_yield_data            = np.loadtxt('output/yield_phi_d_recon_data_2H_exc_Wphi.txt')[:,6]
# phi_d_2H_Wphi_yield_sim             = np.loadtxt('output/yield_phi_d_recon_sim_2H_exc_Wphi.txt')[:,6]
# phi_d_2H_Wphi_yield_tagged          = np.loadtxt('output/yield_phi_d_thrown_tagged_2H_Wphi.txt')[:,6]
# phi_d_2H_Wphi_yield_data_statserr   = np.sqrt(phi_d_2H_Wphi_yield_data)
# phi_d_2H_Wphi_yield_sim_statserr    = np.sqrt(phi_d_2H_Wphi_yield_sim)
# phi_d_2H_Wphi_yield_tagged_statserr = np.sqrt(phi_d_2H_Wphi_yield_tagged)

# # Calculate the efficiency
# phi_d_2H_Wphi_efficiency            = phi_d_2H_Wphi_yield_sim/phi_d_2H_Wphi_yield_tagged
# phi_d_2H_Wphi_efficiency_statserr   = phi_d_2H_Wphi_efficiency * np.sqrt(1/phi_d_2H_Wphi_yield_sim + 1/phi_d_2H_Wphi_yield_tagged)

# # Calculate the results
# phi_d_2H_Wphi_results               = phi_d_2H_Wphi_yield_data/phi_d_2H_Wphi_efficiency  # raw results
# phi_d_2H_Wphi_results               = normalize_distribution(phi_d_2H_Wphi_results, phi_d_2H_Wphi_energy_low, phi_d_2H_Wphi_minust_low) # normalize to have the sum equal to 1
# phi_d_2H_Wphi_results               = 2*np.pi*phi_d_2H_Wphi_results/((phi_d_2H_Wphi_phi_high - phi_d_2H_Wphi_phi_low)/180*np.pi)  # normalize to have the integral equal to 2pi
# phi_d_2H_Wphi_results_statserr      = phi_d_2H_Wphi_results/np.sqrt(phi_d_2H_Wphi_yield_data)

# # Find the indices for the different energy and t bins
# index = []
# for i in range(len(phi_d_2H_Wphi_results)):
#     if (i == 0):
#         index.append(i)
#     elif (i == len(phi_d_2H_Wphi_results) - 1):
#         index.append(i+1)
#     else:
#         if (phi_d_2H_Wphi_energy_low[i] != phi_d_2H_Wphi_energy_low[i-1]) or (phi_d_2H_Wphi_minust_low[i] != phi_d_2H_Wphi_minust_low[i-1]):
#             index.append(i)

# # Plot the data yield
# fig = plt.figure(figsize=(6*(len(index) - 1), 6))
# gs = fig.add_gridspec(1, len(index) - 1, wspace=0)
# axs = gs.subplots(sharex=True, sharey=True)
# for i in range(len(index) - 1):
#     axs[i].errorbar(phi_d_2H_Wphi_phi_center[index[i]:index[i+1]], phi_d_2H_Wphi_yield_data[index[i]:index[i+1]], xerr=phi_d_2H_Wphi_phi_width[index[i]:index[i+1]], yerr=phi_d_2H_Wphi_yield_data_statserr[index[i]:index[i+1]], fmt='k.', label='This work')
#     axs[i].set_title(r'$%.1f<\mathrm{E}_{\gamma}<%.1f\ \mathrm{GeV},\  %.1f<-t<%.1f\ \mathrm{GeV}^2$' % (phi_d_2H_Wphi_energy_low[index[i]], phi_d_2H_Wphi_energy_high[index[i]], phi_d_2H_Wphi_minust_low[index[i]], phi_d_2H_Wphi_minust_high[index[i]]))
# axs[0].set_xlim(-180, 180)
# axs[0].set_ylim(0, 150)
# axs[0].set_xticks(np.arange(-120, 180, 60))
# axs[0].set_ylabel(r'$\mathrm{Yield}$')
# fig.suptitle(r"$d(\gamma, \phi d')$ yield vs $\varphi$")
# fig.supxlabel(r'$\varphi\ [\mathrm{deg}]$')
# plt.legend()
# plt.savefig('output/fig_phi_d_2H_Wphi_yield.png', dpi=300)
# plt.close()

# # Plot the efficiency
# fig = plt.figure(figsize=(6*(len(index) - 1), 6))
# gs = fig.add_gridspec(1, len(index) - 1, wspace=0)
# axs = gs.subplots(sharex=True, sharey=True)
# for i in range(len(index) - 1):
#     axs[i].errorbar(phi_d_2H_Wphi_phi_center[index[i]:index[i+1]], phi_d_2H_Wphi_efficiency[index[i]:index[i+1]], xerr=phi_d_2H_Wphi_phi_width[index[i]:index[i+1]], yerr=phi_d_2H_Wphi_efficiency_statserr[index[i]:index[i+1]], fmt='k.', label='This work')
#     axs[i].set_title(r'$%.1f<\mathrm{E}_{\gamma}<%.1f\ \mathrm{GeV},\  %.1f<-t<%.1f\ \mathrm{GeV}^2$' % (phi_d_2H_Wphi_energy_low[index[i]], phi_d_2H_Wphi_energy_high[index[i]], phi_d_2H_Wphi_minust_low[index[i]], phi_d_2H_Wphi_minust_high[index[i]]))
# axs[0].set_xlim(-180, 180)
# axs[0].set_ylim(0, 0.1)
# axs[0].set_xticks(np.arange(-120, 180, 60))
# axs[0].set_ylabel(r'$\mathrm{Efficiency}$')
# fig.suptitle(r"$d(\gamma, \phi d')$ efficiency vs $\varphi$")
# fig.supxlabel(r'$\varphi\ [\mathrm{deg}]$')
# plt.legend()
# plt.savefig('output/fig_phi_d_2H_Wphi_efficiency.png', dpi=300)
# plt.close()

# # Plot the results
# fig = plt.figure(figsize=(6*(len(index) - 1), 6))
# gs = fig.add_gridspec(1, len(index) - 1, wspace=0)
# axs = gs.subplots(sharex=True, sharey=True)
# for i in range(len(index) - 1):
#     axs[i].errorbar(phi_d_2H_Wphi_phi_center[index[i]:index[i+1]], phi_d_2H_Wphi_results[index[i]:index[i+1]], xerr=phi_d_2H_Wphi_phi_width[index[i]:index[i+1]], yerr=phi_d_2H_Wphi_results_statserr[index[i]:index[i+1]], fmt='k.', label='This work')
#     axs[i].set_title(r'$%.1f<\mathrm{E}_{\gamma}<%.1f\ \mathrm{GeV},\  %.1f<-t<%.1f\ \mathrm{GeV}^2$' % (phi_d_2H_Wphi_energy_low[index[i]], phi_d_2H_Wphi_energy_high[index[i]], phi_d_2H_Wphi_minust_low[index[i]], phi_d_2H_Wphi_minust_high[index[i]]))
#     curve_fit_params, curve_fit_cov = curve_fit(Wphi_func, phi_d_2H_Wphi_phi_center[index[i]:index[i+1]], phi_d_2H_Wphi_results[index[i]:index[i+1]], p0=[0.0])
#     curve_fit_residuals = phi_d_2H_Wphi_results[index[i]:index[i+1]] - Wphi_func(phi_d_2H_Wphi_phi_center[index[i]:index[i+1]], curve_fit_params[0])
#     reduced_chi2 = np.sum((curve_fit_residuals/phi_d_2H_Wphi_results_statserr[index[i]:index[i+1]])**2)/(len(phi_d_2H_Wphi_results[index[i]:index[i+1]])-1)
#     axs[i].plot(np.linspace(-180, 180, 360), Wphi_func(np.linspace(-180, 180, 360), curve_fit_params[0]), 'b--', label='Fit')
#     axs[i].text(-175, 1.95, r'$\rho^0_{1-1}=%.2f\pm%.2f$' % (curve_fit_params[0], np.sqrt(curve_fit_cov[0])), fontsize=10, color='b', ha='left', va='top')
#     axs[i].text(-175, 1.85, r'$\chi^2/\mathrm{dof}=%.2f$' % reduced_chi2, fontsize=10, color='b', ha='left', va='top')
#     axs[i].plot(np.linspace(-180, 180, 360), np.ones(360), 'r--', label='SCHC+NPE')
# axs[0].set_xlim(-180, 180)
# axs[0].set_ylim(0, 2)
# axs[0].set_xticks(np.arange(-120, 180, 60))
# axs[0].set_ylabel(r'$2\pi W(\varphi)$')
# fig.suptitle(r"$d(\gamma, \phi d')$ normalized distribution of $\varphi$")
# fig.supxlabel(r'$\varphi\ [\mathrm{deg}]$')
# plt.legend()
# plt.savefig('output/fig_phi_d_2H_Wphi_results.png', dpi=300)
# plt.close()

# #======================================================================phi_d_2H_WPhi======================================================================

# # Read the bin edges
# phi_d_2H_WPhi_energy_low            = np.loadtxt('configs/bins_phi_d_2H_WPhi.txt')[:,0]
# phi_d_2H_WPhi_energy_high           = np.loadtxt('configs/bins_phi_d_2H_WPhi.txt')[:,1]
# phi_d_2H_WPhi_minust_low            = np.loadtxt('configs/bins_phi_d_2H_WPhi.txt')[:,2]
# phi_d_2H_WPhi_minust_high           = np.loadtxt('configs/bins_phi_d_2H_WPhi.txt')[:,3]
# phi_d_2H_WPhi_Phi_low               = np.loadtxt('configs/bins_phi_d_2H_WPhi.txt')[:,4]
# phi_d_2H_WPhi_Phi_high              = np.loadtxt('configs/bins_phi_d_2H_WPhi.txt')[:,5]
# phi_d_2H_WPhi_Phi_center            = (phi_d_2H_WPhi_Phi_high + phi_d_2H_WPhi_Phi_low) / 2
# phi_d_2H_WPhi_Phi_width             = (phi_d_2H_WPhi_Phi_high - phi_d_2H_WPhi_Phi_low) / 2

# # Read the yield numbers
# phi_d_2H_WPhi_yield_data            = np.loadtxt('output/yield_phi_d_recon_data_2H_exc_WPhi.txt')[:,6]
# phi_d_2H_WPhi_yield_sim             = np.loadtxt('output/yield_phi_d_recon_sim_2H_exc_WPhi.txt')[:,6]
# phi_d_2H_WPhi_yield_tagged          = np.loadtxt('output/yield_phi_d_thrown_tagged_2H_WPhi.txt')[:,6]
# phi_d_2H_WPhi_yield_data_statserr   = np.sqrt(phi_d_2H_WPhi_yield_data)
# phi_d_2H_WPhi_yield_sim_statserr    = np.sqrt(phi_d_2H_WPhi_yield_sim)
# phi_d_2H_WPhi_yield_tagged_statserr = np.sqrt(phi_d_2H_WPhi_yield_tagged)

# # Calculate the efficiency
# phi_d_2H_WPhi_efficiency            = phi_d_2H_WPhi_yield_sim/phi_d_2H_WPhi_yield_tagged
# phi_d_2H_WPhi_efficiency_statserr   = phi_d_2H_WPhi_efficiency * np.sqrt(1/phi_d_2H_WPhi_yield_sim + 1/phi_d_2H_WPhi_yield_tagged)

# # Calculate the results
# phi_d_2H_WPhi_results               = phi_d_2H_WPhi_yield_data/phi_d_2H_WPhi_efficiency  # raw results
# phi_d_2H_WPhi_results               = normalize_distribution(phi_d_2H_WPhi_results, phi_d_2H_WPhi_energy_low, phi_d_2H_WPhi_minust_low) # normalize to have the sum equal to 1
# phi_d_2H_WPhi_results               = 2*np.pi*phi_d_2H_WPhi_results/((phi_d_2H_WPhi_Phi_high - phi_d_2H_WPhi_Phi_low)/180*np.pi)  # normalize to have the integral equal to 2pi
# phi_d_2H_WPhi_results_statserr      = phi_d_2H_WPhi_results/np.sqrt(phi_d_2H_WPhi_yield_data)

# # Find the indices for the different energy and t bins
# index = []
# for i in range(len(phi_d_2H_WPhi_results)):
#     if (i == 0):
#         index.append(i)
#     elif (i == len(phi_d_2H_WPhi_results) - 1):
#         index.append(i+1)
#     else:
#         if (phi_d_2H_WPhi_energy_low[i] != phi_d_2H_WPhi_energy_low[i-1]) or (phi_d_2H_WPhi_minust_low[i] != phi_d_2H_WPhi_minust_low[i-1]):
#             index.append(i)

# # Plot the data yield
# fig = plt.figure(figsize=(6*(len(index) - 1), 6))
# gs = fig.add_gridspec(1, len(index) - 1, wspace=0)
# axs = gs.subplots(sharex=True, sharey=True)
# for i in range(len(index) - 1):
#     axs[i].errorbar(phi_d_2H_WPhi_Phi_center[index[i]:index[i+1]], phi_d_2H_WPhi_yield_data[index[i]:index[i+1]], xerr=phi_d_2H_WPhi_Phi_width[index[i]:index[i+1]], yerr=phi_d_2H_WPhi_yield_data_statserr[index[i]:index[i+1]], fmt='k.', label='This work')
#     axs[i].set_title(r'$%.1f<\mathrm{E}_{\gamma}<%.1f\ \mathrm{GeV},\  %.1f<-t<%.1f\ \mathrm{GeV}^2$' % (phi_d_2H_WPhi_energy_low[index[i]], phi_d_2H_WPhi_energy_high[index[i]], phi_d_2H_WPhi_minust_low[index[i]], phi_d_2H_WPhi_minust_high[index[i]]))
# axs[0].set_xlim(-180, 180)
# axs[0].set_ylim(0, 150)
# axs[0].set_xticks(np.arange(-120, 180, 60))
# axs[0].set_ylabel(r'$\mathrm{Yield}$')
# fig.suptitle(r"$d(\gamma, \phi d')$ yield vs $\Phi$")
# fig.supxlabel(r'$\Phi\ [\mathrm{deg}]$')
# plt.legend()
# plt.savefig('output/fig_phi_d_2H_WPhi_yield.png', dpi=300)
# plt.close()

# # Plot the efficiency
# fig = plt.figure(figsize=(6*(len(index) - 1), 6))
# gs = fig.add_gridspec(1, len(index) - 1, wspace=0)
# axs = gs.subplots(sharex=True, sharey=True)
# for i in range(len(index) - 1):
#     axs[i].errorbar(phi_d_2H_WPhi_Phi_center[index[i]:index[i+1]], phi_d_2H_WPhi_efficiency[index[i]:index[i+1]], xerr=phi_d_2H_WPhi_Phi_width[index[i]:index[i+1]], yerr=phi_d_2H_WPhi_efficiency_statserr[index[i]:index[i+1]], fmt='k.', label='This work')
#     axs[i].set_title(r'$%.1f<\mathrm{E}_{\gamma}<%.1f\ \mathrm{GeV},\  %.1f<-t<%.1f\ \mathrm{GeV}^2$' % (phi_d_2H_WPhi_energy_low[index[i]], phi_d_2H_WPhi_energy_high[index[i]], phi_d_2H_WPhi_minust_low[index[i]], phi_d_2H_WPhi_minust_high[index[i]]))
# axs[0].set_xlim(-180, 180)
# axs[0].set_ylim(0, 0.1)
# axs[0].set_xticks(np.arange(-120, 180, 60))
# axs[0].set_ylabel(r'$\mathrm{Efficiency}$')
# fig.suptitle(r"$d(\gamma, \phi d')$ efficiency vs $\Phi$")
# fig.supxlabel(r'$\Phi\ [\mathrm{deg}]$')
# plt.legend()
# plt.savefig('output/fig_phi_d_2H_WPhi_efficiency.png', dpi=300)
# plt.close()

# # Plot the results
# fig = plt.figure(figsize=(6*(len(index) - 1), 6))
# gs = fig.add_gridspec(1, len(index) - 1, wspace=0)
# axs = gs.subplots(sharex=True, sharey=True)
# for i in range(len(index) - 1):
#     axs[i].errorbar(phi_d_2H_WPhi_Phi_center[index[i]:index[i+1]], phi_d_2H_WPhi_results[index[i]:index[i+1]], xerr=phi_d_2H_WPhi_Phi_width[index[i]:index[i+1]], yerr=phi_d_2H_WPhi_results_statserr[index[i]:index[i+1]], fmt='k.', label='This work')
#     axs[i].set_title(r'$%.1f<\mathrm{E}_{\gamma}<%.1f\ \mathrm{GeV},\  %.1f<-t<%.1f\ \mathrm{GeV}^2$' % (phi_d_2H_WPhi_energy_low[index[i]], phi_d_2H_WPhi_energy_high[index[i]], phi_d_2H_WPhi_minust_low[index[i]], phi_d_2H_WPhi_minust_high[index[i]]))
#     curve_fit_params, curve_fit_cov = curve_fit(Wphi_func, phi_d_2H_WPhi_Phi_center[index[i]:index[i+1]], phi_d_2H_WPhi_results[index[i]:index[i+1]], p0=[0.0])
#     curve_fit_residuals = phi_d_2H_WPhi_results[index[i]:index[i+1]] - Wphi_func(phi_d_2H_WPhi_Phi_center[index[i]:index[i+1]], curve_fit_params[0])
#     reduced_chi2 = np.sum((curve_fit_residuals/phi_d_2H_WPhi_results_statserr[index[i]:index[i+1]])**2)/(len(phi_d_2H_WPhi_results[index[i]:index[i+1]])-1)
#     axs[i].plot(np.linspace(-180, 180, 360), Wphi_func(np.linspace(-180, 180, 360), curve_fit_params[0]), 'b--', label='Fit')
#     axs[i].text(-175, 1.95, r'$\rho^1_{00}=%.2f\pm%.2f$' % (2*curve_fit_params[0]/3/0.3, np.sqrt(curve_fit_cov[0])), fontsize=10, color='b', ha='left', va='top')
#     axs[i].text(-175, 1.85, r'$\chi^2/\mathrm{dof}=%.2f$' % reduced_chi2, fontsize=10, color='b', ha='left', va='top')
#     axs[i].plot(np.linspace(-180, 180, 360), np.ones(360), 'r--', label='SCHC+NPE')
# axs[0].set_xlim(-180, 180)
# axs[0].set_ylim(0, 2)
# axs[0].set_xticks(np.arange(-120, 180, 60))
# axs[0].set_ylabel(r'$2\pi W(\Phi)$')
# fig.suptitle(r"$d(\gamma, \phi d')$ normalized distribution of $\Phi$")
# fig.supxlabel(r'$\Phi\ [\mathrm{deg}]$')
# plt.legend()
# plt.savefig('output/fig_phi_d_2H_WPhi_results.png', dpi=300)
# plt.close()

# #======================================================================phi_d_2H_Wpsi======================================================================

# # Read the bin edges
# phi_d_2H_Wpsi_energy_low            = np.loadtxt('configs/bins_phi_d_2H_Wpsi.txt')[:,0]
# phi_d_2H_Wpsi_energy_high           = np.loadtxt('configs/bins_phi_d_2H_Wpsi.txt')[:,1]
# phi_d_2H_Wpsi_minust_low            = np.loadtxt('configs/bins_phi_d_2H_Wpsi.txt')[:,2]
# phi_d_2H_Wpsi_minust_high           = np.loadtxt('configs/bins_phi_d_2H_Wpsi.txt')[:,3]
# phi_d_2H_Wpsi_psi_low               = np.loadtxt('configs/bins_phi_d_2H_Wpsi.txt')[:,4]
# phi_d_2H_Wpsi_psi_high              = np.loadtxt('configs/bins_phi_d_2H_Wpsi.txt')[:,5]
# phi_d_2H_Wpsi_psi_center            = (phi_d_2H_Wpsi_psi_high + phi_d_2H_Wpsi_psi_low) / 2
# phi_d_2H_Wpsi_psi_width             = (phi_d_2H_Wpsi_psi_high - phi_d_2H_Wpsi_psi_low) / 2

# # Read the yield numbers
# phi_d_2H_Wpsi_yield_data            = np.loadtxt('output/yield_phi_d_recon_data_2H_exc_Wpsi.txt')[:,6]
# phi_d_2H_Wpsi_yield_sim             = np.loadtxt('output/yield_phi_d_recon_sim_2H_exc_Wpsi.txt')[:,6]
# phi_d_2H_Wpsi_yield_tagged          = np.loadtxt('output/yield_phi_d_thrown_tagged_2H_Wpsi.txt')[:,6]
# phi_d_2H_Wpsi_yield_data_statserr   = np.sqrt(phi_d_2H_Wpsi_yield_data)
# phi_d_2H_Wpsi_yield_sim_statserr    = np.sqrt(phi_d_2H_Wpsi_yield_sim)
# phi_d_2H_Wpsi_yield_tagged_statserr = np.sqrt(phi_d_2H_Wpsi_yield_tagged)

# # Calculate the efficiency
# phi_d_2H_Wpsi_efficiency            = phi_d_2H_Wpsi_yield_sim/phi_d_2H_Wpsi_yield_tagged
# phi_d_2H_Wpsi_efficiency_statserr   = phi_d_2H_Wpsi_efficiency * np.sqrt(1/phi_d_2H_Wpsi_yield_sim + 1/phi_d_2H_Wpsi_yield_tagged)

# # Calculate the results
# phi_d_2H_Wpsi_results               = phi_d_2H_Wpsi_yield_data/phi_d_2H_Wpsi_efficiency  # raw results
# phi_d_2H_Wpsi_results               = normalize_distribution(phi_d_2H_Wpsi_results, phi_d_2H_Wpsi_energy_low, phi_d_2H_Wpsi_minust_low) # normalize to have the sum equal to 1
# phi_d_2H_Wpsi_results               = 2*np.pi*phi_d_2H_Wpsi_results/((phi_d_2H_Wpsi_psi_high - phi_d_2H_Wpsi_psi_low)/180*np.pi)  # normalize to have the integral equal to 2pi
# phi_d_2H_Wpsi_results_statserr      = phi_d_2H_Wpsi_results/np.sqrt(phi_d_2H_Wpsi_yield_data)

# # Find the indices for the different energy and t bins
# index = []
# for i in range(len(phi_d_2H_Wpsi_results)):
#     if (i == 0):
#         index.append(i)
#     elif (i == len(phi_d_2H_Wpsi_results) - 1):
#         index.append(i+1)
#     else:
#         if (phi_d_2H_Wpsi_energy_low[i] != phi_d_2H_Wpsi_energy_low[i-1]) or (phi_d_2H_Wpsi_minust_low[i] != phi_d_2H_Wpsi_minust_low[i-1]):
#             index.append(i)

# # Plot the data yield
# fig = plt.figure(figsize=(6*(len(index) - 1), 6))
# gs = fig.add_gridspec(1, len(index) - 1, wspace=0)
# axs = gs.subplots(sharex=True, sharey=True)
# for i in range(len(index) - 1):
#     axs[i].errorbar(phi_d_2H_Wpsi_psi_center[index[i]:index[i+1]], phi_d_2H_Wpsi_yield_data[index[i]:index[i+1]], xerr=phi_d_2H_Wpsi_psi_width[index[i]:index[i+1]], yerr=phi_d_2H_Wpsi_yield_data_statserr[index[i]:index[i+1]], fmt='k.', label='This work')
#     axs[i].set_title(r'$%.1f<\mathrm{E}_{\gamma}<%.1f\ \mathrm{GeV},\  %.1f<-t<%.1f\ \mathrm{GeV}^2$' % (phi_d_2H_Wpsi_energy_low[index[i]], phi_d_2H_Wpsi_energy_high[index[i]], phi_d_2H_Wpsi_minust_low[index[i]], phi_d_2H_Wpsi_minust_high[index[i]]))
# axs[0].set_xlim(-180, 180)
# axs[0].set_ylim(0, 150)
# axs[0].set_xticks(np.arange(-120, 180, 60))
# axs[0].set_ylabel(r'$\mathrm{Yield}$')
# fig.suptitle(r"$d(\gamma, \phi d')$ yield vs $\psi$")
# fig.supxlabel(r'$\psi\ [\mathrm{deg}]$')
# plt.legend()
# plt.savefig('output/fig_phi_d_2H_Wpsi_yield.png', dpi=300)
# plt.close()

# # Plot the efficiency
# fig = plt.figure(figsize=(6*(len(index) - 1), 6))
# gs = fig.add_gridspec(1, len(index) - 1, wspace=0)
# axs = gs.subplots(sharex=True, sharey=True)
# for i in range(len(index) - 1):
#     axs[i].errorbar(phi_d_2H_Wpsi_psi_center[index[i]:index[i+1]], phi_d_2H_Wpsi_efficiency[index[i]:index[i+1]], xerr=phi_d_2H_Wpsi_psi_width[index[i]:index[i+1]], yerr=phi_d_2H_Wpsi_efficiency_statserr[index[i]:index[i+1]], fmt='k.', label='This work')
#     axs[i].set_title(r'$%.1f<\mathrm{E}_{\gamma}<%.1f\ \mathrm{GeV},\  %.1f<-t<%.1f\ \mathrm{GeV}^2$' % (phi_d_2H_Wpsi_energy_low[index[i]], phi_d_2H_Wpsi_energy_high[index[i]], phi_d_2H_Wpsi_minust_low[index[i]], phi_d_2H_Wpsi_minust_high[index[i]]))
# axs[0].set_xlim(-180, 180)
# axs[0].set_ylim(0, 0.1)
# axs[0].set_xticks(np.arange(-120, 180, 60))
# axs[0].set_ylabel(r'$\mathrm{Efficiency}$')
# fig.suptitle(r"$d(\gamma, \phi d')$ efficiency vs $\psi$")
# fig.supxlabel(r'$\psi\ [\mathrm{deg}]$')
# plt.legend()
# plt.savefig('output/fig_phi_d_2H_Wpsi_efficiency.png', dpi=300)
# plt.close()

# # Plot the results
# fig = plt.figure(figsize=(6*(len(index) - 1), 6))
# gs = fig.add_gridspec(1, len(index) - 1, wspace=0)
# axs = gs.subplots(sharex=True, sharey=True)
# for i in range(len(index) - 1):
#     axs[i].errorbar(phi_d_2H_Wpsi_psi_center[index[i]:index[i+1]], phi_d_2H_Wpsi_results[index[i]:index[i+1]], xerr=phi_d_2H_Wpsi_psi_width[index[i]:index[i+1]], yerr=phi_d_2H_Wpsi_results_statserr[index[i]:index[i+1]], fmt='k.', label='This work')
#     axs[i].set_title(r'$%.1f<\mathrm{E}_{\gamma}<%.1f\ \mathrm{GeV},\  %.1f<-t<%.1f\ \mathrm{GeV}^2$' % (phi_d_2H_Wpsi_energy_low[index[i]], phi_d_2H_Wpsi_energy_high[index[i]], phi_d_2H_Wpsi_minust_low[index[i]], phi_d_2H_Wpsi_minust_high[index[i]]))
#     curve_fit_params, curve_fit_cov = curve_fit(Wphi_func, phi_d_2H_Wpsi_psi_center[index[i]:index[i+1]], phi_d_2H_Wpsi_results[index[i]:index[i+1]], p0=[0.0])
#     curve_fit_residuals = phi_d_2H_Wpsi_results[index[i]:index[i+1]] - Wphi_func(phi_d_2H_Wpsi_psi_center[index[i]:index[i+1]], curve_fit_params[0])
#     reduced_chi2 = np.sum((curve_fit_residuals/phi_d_2H_Wpsi_results_statserr[index[i]:index[i+1]])**2)/(len(phi_d_2H_Wpsi_results[index[i]:index[i+1]])-1)
#     axs[i].plot(np.linspace(-180, 180, 360), Wphi_func(np.linspace(-180, 180, 360), curve_fit_params[0]), 'b--', label='Fit')
#     axs[i].text(-175, 1.95, r'$\rho^1_{1-1}=%.2f\pm%.2f$' % (-curve_fit_params[0]/0.3, np.sqrt(curve_fit_cov[0])), fontsize=10, color='b', ha='left', va='top')
#     axs[i].text(-175, 1.85, r'$\chi^2/\mathrm{dof}=%.2f$' % reduced_chi2, fontsize=10, color='b', ha='left', va='top')
#     axs[i].plot(np.linspace(-180, 180, 360), np.ones(360)+0.3*np.cos(2*np.linspace(-180, 180, 360)/rad_to_deg), 'r--', label='SCHC+NPE')
# axs[0].set_xlim(-180, 180)
# axs[0].set_ylim(0, 2)
# axs[0].set_xticks(np.arange(-120, 180, 60))
# axs[0].set_ylabel(r'$2\pi W(\psi)$')
# fig.suptitle(r"$d(\gamma, \phi d')$ normalized distribution of $\psi$")
# fig.supxlabel(r'$\psi\ [\mathrm{deg}]$')
# plt.legend()
# plt.savefig('output/fig_phi_d_2H_Wpsi_results.png', dpi=300)
# plt.close()

# #======================================================================rho_d_2H_dsdt======================================================================

# # Read the bin edges
# rho_d_2H_dsdt_energy_low            = np.loadtxt('configs/bins_rho_d_2H_dsdt.txt')[:,0]
# rho_d_2H_dsdt_energy_high           = np.loadtxt('configs/bins_rho_d_2H_dsdt.txt')[:,1]
# rho_d_2H_dsdt_minust_low            = np.loadtxt('configs/bins_rho_d_2H_dsdt.txt')[:,2]
# rho_d_2H_dsdt_minust_high           = np.loadtxt('configs/bins_rho_d_2H_dsdt.txt')[:,3]
# rho_d_2H_dsdt_minust_center         = (rho_d_2H_dsdt_minust_high + rho_d_2H_dsdt_minust_low) / 2
# rho_d_2H_dsdt_minust_width          = (rho_d_2H_dsdt_minust_high - rho_d_2H_dsdt_minust_low) / 2

# # Read the yield numbers
# rho_d_2H_dsdt_yield_data            = np.loadtxt('output/yield_rho_d_recon_data_2H_exc_dsdt.txt')[:,4]
# rho_d_2H_dsdt_yield_sim             = np.loadtxt('output/yield_rho_d_recon_sim_2H_exc_dsdt.txt')[:,4]
# rho_d_2H_dsdt_yield_tagged          = np.loadtxt('output/yield_rho_d_thrown_tagged_2H_dsdt.txt')[:,4]
# rho_d_2H_dsdt_yield_data_statserr   = np.sqrt(rho_d_2H_dsdt_yield_data)
# rho_d_2H_dsdt_yield_sim_statserr    = np.sqrt(rho_d_2H_dsdt_yield_sim)
# rho_d_2H_dsdt_yield_tagged_statserr = np.sqrt(rho_d_2H_dsdt_yield_tagged)

# # Calculate the efficiency
# rho_d_2H_dsdt_efficiency            = rho_d_2H_dsdt_yield_sim/rho_d_2H_dsdt_yield_tagged
# rho_d_2H_dsdt_efficiency_statserr   = rho_d_2H_dsdt_efficiency * np.sqrt(1/rho_d_2H_dsdt_yield_sim + 1/rho_d_2H_dsdt_yield_tagged)

# # Calculate the results
# rho_d_2H_dsdt_results               = rho_d_2H_dsdt_yield_data/rho_d_2H_dsdt_efficiency/lumi(rho_d_2H_dsdt_energy_low, rho_d_2H_dsdt_energy_high, '2H')/(rho_d_2H_dsdt_minust_high-rho_d_2H_dsdt_minust_low)/0.489/1000
# rho_d_2H_dsdt_results_statserr      = rho_d_2H_dsdt_results/np.sqrt(rho_d_2H_dsdt_yield_data)

# # Find the indices for the different energy and t bins
# index = []
# for i in range(len(rho_d_2H_dsdt_results)):
#     if (i == 0):
#         index.append(i)
#     elif (i == len(rho_d_2H_dsdt_results) - 1):
#         index.append(i+1)
#     else:
#         if (rho_d_2H_dsdt_energy_low[i] != rho_d_2H_dsdt_energy_low[i-1]):
#             index.append(i)

# # Plot the data yield
# fig = plt.figure(figsize=(8, 6))
# for i in range(len(index) - 1):
#     plt.errorbar(rho_d_2H_dsdt_minust_center[index[i]:index[i+1]], rho_d_2H_dsdt_yield_data[index[i]:index[i+1]], xerr=rho_d_2H_dsdt_minust_width[index[i]:index[i+1]], yerr=rho_d_2H_dsdt_yield_data_statserr[index[i]:index[i+1]], fmt='.', label=r'$%.1f<\mathrm{E}_{\gamma}<%.1f\ \mathrm{GeV}$' % (rho_d_2H_dsdt_energy_low[index[i]], rho_d_2H_dsdt_energy_high[index[i]]))
# plt.title(r"$d(\gamma, \rho d')$ yield vs $-t$")
# plt.xlabel(r'$-t[GeV^2/c]$')
# plt.ylabel(r'$\mathrm{Yield}$')
# plt.xlim(0, 2)
# plt.ylim(0, 300)
# plt.legend()
# plt.savefig('output/fig_rho_d_2H_dsdt_yield.png', dpi=300)
# plt.close()

# # Plot the efficiency
# fig = plt.figure(figsize=(8, 6))
# for i in range(len(index) - 1):
#     plt.errorbar(rho_d_2H_dsdt_minust_center[index[i]:index[i+1]], rho_d_2H_dsdt_efficiency[index[i]:index[i+1]], xerr=rho_d_2H_dsdt_minust_width[index[i]:index[i+1]], yerr=rho_d_2H_dsdt_efficiency_statserr[index[i]:index[i+1]], fmt='.', label=r'$%.1f<\mathrm{E}_{\gamma}<%.1f\ \mathrm{GeV}$' % (rho_d_2H_dsdt_energy_low[index[i]], rho_d_2H_dsdt_energy_high[index[i]]))
# plt.title(r"$d(\gamma, \rho d')$ efficiency vs $-t$")
# plt.xlabel(r'$-t[GeV^2/c]$')
# plt.ylabel(r'$\mathrm{Efficiency}$')
# plt.xlim(0, 2)
# plt.ylim(0, 0.1)
# plt.legend()
# plt.savefig('output/fig_rho_d_2H_dsdt_efficiency.png', dpi=300)
# plt.close()

# # Plot the results
# fig = plt.figure(figsize=(8, 6))
# for i in range(len(index) - 1):
#     plt.errorbar(rho_d_2H_dsdt_minust_center[index[i]:index[i+1]], rho_d_2H_dsdt_results[index[i]:index[i+1]], xerr=rho_d_2H_dsdt_minust_width[index[i]:index[i+1]], yerr=rho_d_2H_dsdt_results_statserr[index[i]:index[i+1]], fmt='.', label=r'$%.1f<\mathrm{E}_{\gamma}<%.1f\ \mathrm{GeV}$' % (rho_d_2H_dsdt_energy_low[index[i]], rho_d_2H_dsdt_energy_high[index[i]]))
# plt.text(0.3, 0.15, 'preliminary', fontsize=15, color='r', style='italic', ha='center', va='center')
# plt.title(r"$d(\gamma, \rho d')$ differential cross section vs $-t$")
# plt.xlabel(r'$-t\ [GeV^2/c]$')
# plt.ylabel(r'$d\sigma/dt\ [nb/(GeV^2/c)]$')
# plt.xlim(0, 2)
# plt.ylim(1e0, 1e4)
# plt.yscale('log')
# plt.legend()
# plt.savefig('output/fig_rho_d_2H_dsdt_results.png', dpi=300)
# plt.close()

# #======================================================================rho_d_2H_Wcostheta======================================================================

# # Read the bin edges
# rho_d_2H_Wcostheta_energy_low               = np.loadtxt('configs/bins_rho_d_2H_Wcostheta.txt')[:,0]
# rho_d_2H_Wcostheta_energy_high              = np.loadtxt('configs/bins_rho_d_2H_Wcostheta.txt')[:,1]
# rho_d_2H_Wcostheta_minust_low               = np.loadtxt('configs/bins_rho_d_2H_Wcostheta.txt')[:,2]
# rho_d_2H_Wcostheta_minust_high              = np.loadtxt('configs/bins_rho_d_2H_Wcostheta.txt')[:,3]
# rho_d_2H_Wcostheta_costheta_low             = np.loadtxt('configs/bins_rho_d_2H_Wcostheta.txt')[:,4]
# rho_d_2H_Wcostheta_costheta_high            = np.loadtxt('configs/bins_rho_d_2H_Wcostheta.txt')[:,5]
# rho_d_2H_Wcostheta_costheta_center          = (rho_d_2H_Wcostheta_costheta_high + rho_d_2H_Wcostheta_costheta_low) / 2
# rho_d_2H_Wcostheta_costheta_width           = (rho_d_2H_Wcostheta_costheta_high - rho_d_2H_Wcostheta_costheta_low) / 2

# # Read the yield numbers
# rho_d_2H_Wcostheta_yield_data               = np.loadtxt('output/yield_rho_d_recon_data_2H_exc_Wcostheta.txt')[:,6]
# rho_d_2H_Wcostheta_yield_sim                = np.loadtxt('output/yield_rho_d_recon_sim_2H_exc_Wcostheta.txt')[:,6]
# rho_d_2H_Wcostheta_yield_tagged             = np.loadtxt('output/yield_rho_d_thrown_tagged_2H_Wcostheta.txt')[:,6]
# rho_d_2H_Wcostheta_yield_data_statserr      = np.sqrt(rho_d_2H_Wcostheta_yield_data)
# rho_d_2H_Wcostheta_yield_sim_statserr       = np.sqrt(rho_d_2H_Wcostheta_yield_sim)
# rho_d_2H_Wcostheta_yield_tagged_statserr    = np.sqrt(rho_d_2H_Wcostheta_yield_tagged)

# # Calculate the efficiency
# rho_d_2H_Wcostheta_efficiency               = rho_d_2H_Wcostheta_yield_sim/rho_d_2H_Wcostheta_yield_tagged
# rho_d_2H_Wcostheta_efficiency_statserr      = rho_d_2H_Wcostheta_efficiency * np.sqrt(1/rho_d_2H_Wcostheta_yield_sim + 1/rho_d_2H_Wcostheta_yield_tagged)

# # Calculate the results
# rho_d_2H_Wcostheta_results                  = rho_d_2H_Wcostheta_yield_data/rho_d_2H_Wcostheta_efficiency  # raw results
# rho_d_2H_Wcostheta_results                  = normalize_distribution(rho_d_2H_Wcostheta_results, rho_d_2H_Wcostheta_energy_low, rho_d_2H_Wcostheta_minust_low) # normalize to have the sum equal to 1
# rho_d_2H_Wcostheta_results                  = rho_d_2H_Wcostheta_results/(rho_d_2H_Wcostheta_costheta_high - rho_d_2H_Wcostheta_costheta_low)  # normalize to have the integral equal to 1
# rho_d_2H_Wcostheta_results_statserr         = rho_d_2H_Wcostheta_results/np.sqrt(rho_d_2H_Wcostheta_yield_data)

# # Find the indices for the different energy and t bins
# index = []
# for i in range(len(rho_d_2H_Wcostheta_results)):
#     if (i == 0):
#         index.append(i)
#     elif (i == len(rho_d_2H_Wcostheta_results) - 1):
#         index.append(i+1)
#     else:
#         if (rho_d_2H_Wcostheta_energy_low[i] != rho_d_2H_Wcostheta_energy_low[i-1]) or (rho_d_2H_Wcostheta_minust_low[i] != rho_d_2H_Wcostheta_minust_low[i-1]):
#             index.append(i)
# minust = np.zeros(len(index) - 1, dtype=float)
# minust_err = np.zeros(len(index) - 1, dtype=float)
# sdme = np.zeros(len(index) - 1, dtype=float)
# sdme_err = np.zeros(len(index) - 1, dtype=float)

# # Plot the data yield
# fig = plt.figure(figsize=(6*(len(index) - 1), 6))
# gs = fig.add_gridspec(1, len(index) - 1, wspace=0)
# axs = gs.subplots(sharex=True, sharey=True)
# for i in range(len(index) - 1):
#     axs[i].errorbar(rho_d_2H_Wcostheta_costheta_center[index[i]:index[i+1]], rho_d_2H_Wcostheta_yield_data[index[i]:index[i+1]], xerr=rho_d_2H_Wcostheta_costheta_width[index[i]:index[i+1]], yerr=rho_d_2H_Wcostheta_yield_data_statserr[index[i]:index[i+1]], fmt='k.', label='This work')
#     axs[i].set_title(r'$%.1f<\mathrm{E}_{\gamma}<%.1f\ \mathrm{GeV},\  %.1f<-t<%.1f\ \mathrm{GeV}^2$' % (rho_d_2H_Wcostheta_energy_low[index[i]], rho_d_2H_Wcostheta_energy_high[index[i]], rho_d_2H_Wcostheta_minust_low[index[i]], rho_d_2H_Wcostheta_minust_high[index[i]]))
# axs[0].set_xlim(-1, 1)
# axs[0].set_ylim(0, 150)
# axs[0].set_xticks(np.arange(-0.75, 0.9, 0.25))
# axs[0].set_ylabel(r'$\mathrm{Yield}$')
# fig.suptitle(r"$d(\gamma, \rho d')$ yield vs $\cos\vartheta$")
# fig.supxlabel(r'$\cos\vartheta$')
# plt.legend()
# plt.savefig('output/fig_rho_d_2H_Wcostheta_yield.png', dpi=300)
# plt.close()

# # Plot the efficiency
# fig = plt.figure(figsize=(6*(len(index) - 1), 6))
# gs = fig.add_gridspec(1, len(index) - 1, wspace=0)
# axs = gs.subplots(sharex=True, sharey=True)
# for i in range(len(index) - 1):
#     axs[i].errorbar(rho_d_2H_Wcostheta_costheta_center[index[i]:index[i+1]], rho_d_2H_Wcostheta_efficiency[index[i]:index[i+1]], xerr=rho_d_2H_Wcostheta_costheta_width[index[i]:index[i+1]], yerr=rho_d_2H_Wcostheta_efficiency_statserr[index[i]:index[i+1]], fmt='k.', label='This work')
#     axs[i].set_title(r'$%.1f<\mathrm{E}_{\gamma}<%.1f\ \mathrm{GeV},\  %.1f<-t<%.1f\ \mathrm{GeV}^2$' % (rho_d_2H_Wcostheta_energy_low[index[i]], rho_d_2H_Wcostheta_energy_high[index[i]], rho_d_2H_Wcostheta_minust_low[index[i]], rho_d_2H_Wcostheta_minust_high[index[i]]))
# axs[0].set_xlim(-1, 1)
# axs[0].set_ylim(0, 0.1)
# axs[0].set_xticks(np.arange(-0.75, 0.9, 0.25))
# axs[0].set_ylabel(r'$\mathrm{Efficiency}$')
# fig.suptitle(r"$d(\gamma, \rho d')$ efficiency vs $\cos\vartheta$")
# fig.supxlabel(r'$\cos\vartheta$')
# plt.legend()
# plt.savefig('output/fig_rho_d_2H_Wcostheta_efficiency.png', dpi=300)
# plt.close()

# # Plot the results
# fig = plt.figure(figsize=(6*(len(index) - 1), 6))
# gs = fig.add_gridspec(1, len(index) - 1, wspace=0)
# axs = gs.subplots(sharex=True, sharey=True)
# for i in range(len(index) - 1):
#     axs[i].errorbar(rho_d_2H_Wcostheta_costheta_center[index[i]:index[i+1]], rho_d_2H_Wcostheta_results[index[i]:index[i+1]], xerr=rho_d_2H_Wcostheta_costheta_width[index[i]:index[i+1]], yerr=rho_d_2H_Wcostheta_results_statserr[index[i]:index[i+1]], fmt='k.', label='This work')
#     axs[i].set_title(r'$%.1f<\mathrm{E}_{\gamma}<%.1f\ \mathrm{GeV},\  %.1f<-t<%.1f\ \mathrm{GeV}^2$' % (rho_d_2H_Wcostheta_energy_low[index[i]], rho_d_2H_Wcostheta_energy_high[index[i]], rho_d_2H_Wcostheta_minust_low[index[i]], rho_d_2H_Wcostheta_minust_high[index[i]]))
#     curve_fit_params, curve_fit_cov = curve_fit(Wcostheta_func, rho_d_2H_Wcostheta_costheta_center[index[i]:index[i+1]], rho_d_2H_Wcostheta_results[index[i]:index[i+1]], p0=[0.0, 0.0])
#     curve_fit_residuals = rho_d_2H_Wcostheta_results[index[i]:index[i+1]] - Wcostheta_func(rho_d_2H_Wcostheta_costheta_center[index[i]:index[i+1]], curve_fit_params[0], curve_fit_params[1])
#     reduced_chi2 = np.sum((curve_fit_residuals/rho_d_2H_Wcostheta_results_statserr[index[i]:index[i+1]])**2)/(len(rho_d_2H_Wcostheta_results[index[i]:index[i+1]])-2)
#     minust[i] = (rho_d_2H_Wcostheta_minust_low[index[i]] + rho_d_2H_Wcostheta_minust_high[index[i]]) / 2
#     minust_err[i] = (rho_d_2H_Wcostheta_minust_high[index[i]] - rho_d_2H_Wcostheta_minust_low[index[i]]) / 2
#     sdme[i] = curve_fit_params[0]
#     sdme_err[i] = np.sqrt(curve_fit_cov[0][0])
#     axs[i].plot(np.linspace(-1, 1, 100), Wcostheta_func(np.linspace(-1, 1, 100), curve_fit_params[0], curve_fit_params[1]), 'b--', label='Fit')
#     axs[i].text(-0.9, 0.95, r'$\rho^0_{00}=%.2f\pm%.2f$' % (curve_fit_params[0], np.sqrt(curve_fit_cov[0][0])), fontsize=10, color='b', ha='left', va='top')
#     axs[i].text(-0.9, 0.90, r'$\alpha=%.2e\pm%.2e$' % (curve_fit_params[1], np.sqrt(curve_fit_cov[1][1])), fontsize=10, color='b', ha='left', va='top')
#     axs[i].text(-0.9, 0.85, r'$\chi^2/\mathrm{dof}=%.2f$' % reduced_chi2, fontsize=10, color='b', ha='left', va='top')
#     axs[i].plot(np.linspace(-1, 1, 100), 0.75*(1-np.linspace(-1, 1, 100)**2), 'r--', label='SCHC+NPE')
# axs[0].set_xlim(-1, 1)
# axs[0].set_ylim(0, 1)
# axs[0].set_xticks(np.arange(-0.75, 0.9, 0.25))
# axs[0].set_ylabel(r'$W(\cos\vartheta)$')
# fig.suptitle(r"$d(\gamma, \rho d')$ normalized distribution of $\cos\vartheta$")
# fig.supxlabel(r'$\cos\vartheta$')
# plt.legend()
# plt.savefig('output/fig_rho_d_2H_Wcostheta_results.png', dpi=300)
# plt.close()

# plt.errorbar(minust, sdme, xerr=minust_err, yerr=sdme_err, fmt='k.', label='This work')
# plt.plot(np.linspace(0, 1, 100), np.zeros(100), 'r--')
# plt.xlabel(r'$-t\ [\mathrm{GeV}^2/c^2]$')
# plt.ylabel(r'$\rho^0_{00}$')
# plt.xlim(0, 1)
# plt.ylim(-0.15, 0.15)
# plt.savefig('output/fig_rho_d_2H_Wcostheta_sdme.png', dpi=300)
# plt.close()

# #======================================================================rho_d_2H_Wphi======================================================================

# # Read the bin edges
# rho_d_2H_Wphi_energy_low            = np.loadtxt('configs/bins_rho_d_2H_Wphi.txt')[:,0]
# rho_d_2H_Wphi_energy_high           = np.loadtxt('configs/bins_rho_d_2H_Wphi.txt')[:,1]
# rho_d_2H_Wphi_minust_low            = np.loadtxt('configs/bins_rho_d_2H_Wphi.txt')[:,2]
# rho_d_2H_Wphi_minust_high           = np.loadtxt('configs/bins_rho_d_2H_Wphi.txt')[:,3]
# rho_d_2H_Wphi_phi_low               = np.loadtxt('configs/bins_rho_d_2H_Wphi.txt')[:,4]
# rho_d_2H_Wphi_phi_high              = np.loadtxt('configs/bins_rho_d_2H_Wphi.txt')[:,5]
# rho_d_2H_Wphi_phi_center            = (rho_d_2H_Wphi_phi_high + rho_d_2H_Wphi_phi_low) / 2
# rho_d_2H_Wphi_phi_width             = (rho_d_2H_Wphi_phi_high - rho_d_2H_Wphi_phi_low) / 2

# # Read the yield numbers
# rho_d_2H_Wphi_yield_data            = np.loadtxt('output/yield_rho_d_recon_data_2H_exc_Wphi.txt')[:,6]
# rho_d_2H_Wphi_yield_sim             = np.loadtxt('output/yield_rho_d_recon_sim_2H_exc_Wphi.txt')[:,6]
# rho_d_2H_Wphi_yield_tagged          = np.loadtxt('output/yield_rho_d_thrown_tagged_2H_Wphi.txt')[:,6]
# rho_d_2H_Wphi_yield_data_statserr   = np.sqrt(rho_d_2H_Wphi_yield_data)
# rho_d_2H_Wphi_yield_sim_statserr    = np.sqrt(rho_d_2H_Wphi_yield_sim)
# rho_d_2H_Wphi_yield_tagged_statserr = np.sqrt(rho_d_2H_Wphi_yield_tagged)

# # Calculate the efficiency
# rho_d_2H_Wphi_efficiency            = rho_d_2H_Wphi_yield_sim/rho_d_2H_Wphi_yield_tagged
# rho_d_2H_Wphi_efficiency_statserr   = rho_d_2H_Wphi_efficiency * np.sqrt(1/rho_d_2H_Wphi_yield_sim + 1/rho_d_2H_Wphi_yield_tagged)

# # Calculate the results
# rho_d_2H_Wphi_results               = rho_d_2H_Wphi_yield_data/rho_d_2H_Wphi_efficiency  # raw results
# rho_d_2H_Wphi_results               = normalize_distribution(rho_d_2H_Wphi_results, rho_d_2H_Wphi_energy_low, rho_d_2H_Wphi_minust_low) # normalize to have the sum equal to 1
# rho_d_2H_Wphi_results               = 2*np.pi*rho_d_2H_Wphi_results/((rho_d_2H_Wphi_phi_high - rho_d_2H_Wphi_phi_low)/180*np.pi)  # normalize to have the integral equal to 2pi
# rho_d_2H_Wphi_results_statserr      = rho_d_2H_Wphi_results/np.sqrt(rho_d_2H_Wphi_yield_data)

# # Find the indices for the different energy and t bins
# index = []
# for i in range(len(rho_d_2H_Wphi_results)):
#     if (i == 0):
#         index.append(i)
#     elif (i == len(rho_d_2H_Wphi_results) - 1):
#         index.append(i+1)
#     else:
#         if (rho_d_2H_Wphi_energy_low[i] != rho_d_2H_Wphi_energy_low[i-1]) or (rho_d_2H_Wphi_minust_low[i] != rho_d_2H_Wphi_minust_low[i-1]):
#             index.append(i)

# # Plot the data yield
# fig = plt.figure(figsize=(6*(len(index) - 1), 6))
# gs = fig.add_gridspec(1, len(index) - 1, wspace=0)
# axs = gs.subplots(sharex=True, sharey=True)
# for i in range(len(index) - 1):
#     axs[i].errorbar(rho_d_2H_Wphi_phi_center[index[i]:index[i+1]], rho_d_2H_Wphi_yield_data[index[i]:index[i+1]], xerr=rho_d_2H_Wphi_phi_width[index[i]:index[i+1]], yerr=rho_d_2H_Wphi_yield_data_statserr[index[i]:index[i+1]], fmt='k.', label='This work')
#     axs[i].set_title(r'$%.1f<\mathrm{E}_{\gamma}<%.1f\ \mathrm{GeV},\  %.1f<-t<%.1f\ \mathrm{GeV}^2$' % (rho_d_2H_Wphi_energy_low[index[i]], rho_d_2H_Wphi_energy_high[index[i]], rho_d_2H_Wphi_minust_low[index[i]], rho_d_2H_Wphi_minust_high[index[i]]))
# axs[0].set_xlim(-180, 180)
# axs[0].set_ylim(0, 150)
# axs[0].set_xticks(np.arange(-120, 180, 60))
# axs[0].set_ylabel(r'$\mathrm{Yield}$')
# fig.suptitle(r"$d(\gamma, \rho d')$ yield vs $\varphi$")
# fig.supxlabel(r'$\varphi\ [\mathrm{deg}]$')
# plt.legend()
# plt.savefig('output/fig_rho_d_2H_Wphi_yield.png', dpi=300)
# plt.close()

# # Plot the efficiency
# fig = plt.figure(figsize=(6*(len(index) - 1), 6))
# gs = fig.add_gridspec(1, len(index) - 1, wspace=0)
# axs = gs.subplots(sharex=True, sharey=True)
# for i in range(len(index) - 1):
#     axs[i].errorbar(rho_d_2H_Wphi_phi_center[index[i]:index[i+1]], rho_d_2H_Wphi_efficiency[index[i]:index[i+1]], xerr=rho_d_2H_Wphi_phi_width[index[i]:index[i+1]], yerr=rho_d_2H_Wphi_efficiency_statserr[index[i]:index[i+1]], fmt='k.', label='This work')
#     axs[i].set_title(r'$%.1f<\mathrm{E}_{\gamma}<%.1f\ \mathrm{GeV},\  %.1f<-t<%.1f\ \mathrm{GeV}^2$' % (rho_d_2H_Wphi_energy_low[index[i]], rho_d_2H_Wphi_energy_high[index[i]], rho_d_2H_Wphi_minust_low[index[i]], rho_d_2H_Wphi_minust_high[index[i]]))
# axs[0].set_xlim(-180, 180)
# axs[0].set_ylim(0, 0.1)
# axs[0].set_xticks(np.arange(-120, 180, 60))
# axs[0].set_ylabel(r'$\mathrm{Efficiency}$')
# fig.suptitle(r"$d(\gamma, \rho d')$ efficiency vs $\varphi$")
# fig.supxlabel(r'$\varphi\ [\mathrm{deg}]$')
# plt.legend()
# plt.savefig('output/fig_rho_d_2H_Wphi_efficiency.png', dpi=300)
# plt.close()

# # Plot the results
# fig = plt.figure(figsize=(6*(len(index) - 1), 6))
# gs = fig.add_gridspec(1, len(index) - 1, wspace=0)
# axs = gs.subplots(sharex=True, sharey=True)
# for i in range(len(index) - 1):
#     axs[i].errorbar(rho_d_2H_Wphi_phi_center[index[i]:index[i+1]], rho_d_2H_Wphi_results[index[i]:index[i+1]], xerr=rho_d_2H_Wphi_phi_width[index[i]:index[i+1]], yerr=rho_d_2H_Wphi_results_statserr[index[i]:index[i+1]], fmt='k.', label='This work')
#     axs[i].set_title(r'$%.1f<\mathrm{E}_{\gamma}<%.1f\ \mathrm{GeV},\  %.1f<-t<%.1f\ \mathrm{GeV}^2$' % (rho_d_2H_Wphi_energy_low[index[i]], rho_d_2H_Wphi_energy_high[index[i]], rho_d_2H_Wphi_minust_low[index[i]], rho_d_2H_Wphi_minust_high[index[i]]))
#     curve_fit_params, curve_fit_cov = curve_fit(Wphi_func, rho_d_2H_Wphi_phi_center[index[i]:index[i+1]], rho_d_2H_Wphi_results[index[i]:index[i+1]], p0=[0.0])
#     curve_fit_residuals = rho_d_2H_Wphi_results[index[i]:index[i+1]] - Wphi_func(rho_d_2H_Wphi_phi_center[index[i]:index[i+1]], curve_fit_params[0])
#     reduced_chi2 = np.sum((curve_fit_residuals/rho_d_2H_Wphi_results_statserr[index[i]:index[i+1]])**2)/(len(rho_d_2H_Wphi_results[index[i]:index[i+1]])-1)
#     axs[i].plot(np.linspace(-180, 180, 360), Wphi_func(np.linspace(-180, 180, 360), curve_fit_params[0]), 'b--', label='Fit')
#     axs[i].text(-175, 1.95, r'$\rho^0_{1-1}=%.2f\pm%.2f$' % (curve_fit_params[0], np.sqrt(curve_fit_cov[0])), fontsize=10, color='b', ha='left', va='top')
#     axs[i].text(-175, 1.85, r'$\chi^2/\mathrm{dof}=%.2f$' % reduced_chi2, fontsize=10, color='b', ha='left', va='top')
#     axs[i].plot(np.linspace(-180, 180, 360), np.ones(360), 'r--', label='SCHC+NPE')
# axs[0].set_xlim(-180, 180)
# axs[0].set_ylim(0, 2)
# axs[0].set_xticks(np.arange(-120, 180, 60))
# axs[0].set_ylabel(r'$2\pi W(\varphi)$')
# fig.suptitle(r"$d(\gamma, \rho d')$ normalized distribution of $\varphi$")
# fig.supxlabel(r'$\varphi\ [\mathrm{deg}]$')
# plt.legend()
# plt.savefig('output/fig_rho_d_2H_Wphi_results.png', dpi=300)
# plt.close()

# #======================================================================rho_d_2H_WPhi======================================================================

# # Read the bin edges
# rho_d_2H_WPhi_energy_low            = np.loadtxt('configs/bins_rho_d_2H_WPhi.txt')[:,0]
# rho_d_2H_WPhi_energy_high           = np.loadtxt('configs/bins_rho_d_2H_WPhi.txt')[:,1]
# rho_d_2H_WPhi_minust_low            = np.loadtxt('configs/bins_rho_d_2H_WPhi.txt')[:,2]
# rho_d_2H_WPhi_minust_high           = np.loadtxt('configs/bins_rho_d_2H_WPhi.txt')[:,3]
# rho_d_2H_WPhi_Phi_low               = np.loadtxt('configs/bins_rho_d_2H_WPhi.txt')[:,4]
# rho_d_2H_WPhi_Phi_high              = np.loadtxt('configs/bins_rho_d_2H_WPhi.txt')[:,5]
# rho_d_2H_WPhi_Phi_center            = (rho_d_2H_WPhi_Phi_high + rho_d_2H_WPhi_Phi_low) / 2
# rho_d_2H_WPhi_Phi_width             = (rho_d_2H_WPhi_Phi_high - rho_d_2H_WPhi_Phi_low) / 2

# # Read the yield numbers
# rho_d_2H_WPhi_yield_data            = np.loadtxt('output/yield_rho_d_recon_data_2H_exc_WPhi.txt')[:,6]
# rho_d_2H_WPhi_yield_sim             = np.loadtxt('output/yield_rho_d_recon_sim_2H_exc_WPhi.txt')[:,6]
# rho_d_2H_WPhi_yield_tagged          = np.loadtxt('output/yield_rho_d_thrown_tagged_2H_WPhi.txt')[:,6]
# rho_d_2H_WPhi_yield_data_statserr   = np.sqrt(rho_d_2H_WPhi_yield_data)
# rho_d_2H_WPhi_yield_sim_statserr    = np.sqrt(rho_d_2H_WPhi_yield_sim)
# rho_d_2H_WPhi_yield_tagged_statserr = np.sqrt(rho_d_2H_WPhi_yield_tagged)

# # Calculate the efficiency
# rho_d_2H_WPhi_efficiency            = rho_d_2H_WPhi_yield_sim/rho_d_2H_WPhi_yield_tagged
# rho_d_2H_WPhi_efficiency_statserr   = rho_d_2H_WPhi_efficiency * np.sqrt(1/rho_d_2H_WPhi_yield_sim + 1/rho_d_2H_WPhi_yield_tagged)

# # Calculate the results
# rho_d_2H_WPhi_results               = rho_d_2H_WPhi_yield_data/rho_d_2H_WPhi_efficiency  # raw results
# rho_d_2H_WPhi_results               = normalize_distribution(rho_d_2H_WPhi_results, rho_d_2H_WPhi_energy_low, rho_d_2H_WPhi_minust_low) # normalize to have the sum equal to 1
# rho_d_2H_WPhi_results               = 2*np.pi*rho_d_2H_WPhi_results/((rho_d_2H_WPhi_Phi_high - rho_d_2H_WPhi_Phi_low)/180*np.pi)  # normalize to have the integral equal to 2pi
# rho_d_2H_WPhi_results_statserr      = rho_d_2H_WPhi_results/np.sqrt(rho_d_2H_WPhi_yield_data)

# # Find the indices for the different energy and t bins
# index = []
# for i in range(len(rho_d_2H_WPhi_results)):
#     if (i == 0):
#         index.append(i)
#     elif (i == len(rho_d_2H_WPhi_results) - 1):
#         index.append(i+1)
#     else:
#         if (rho_d_2H_WPhi_energy_low[i] != rho_d_2H_WPhi_energy_low[i-1]) or (rho_d_2H_WPhi_minust_low[i] != rho_d_2H_WPhi_minust_low[i-1]):
#             index.append(i)

# # Plot the data yield
# fig = plt.figure(figsize=(6*(len(index) - 1), 6))
# gs = fig.add_gridspec(1, len(index) - 1, wspace=0)
# axs = gs.subplots(sharex=True, sharey=True)
# for i in range(len(index) - 1):
#     axs[i].errorbar(rho_d_2H_WPhi_Phi_center[index[i]:index[i+1]], rho_d_2H_WPhi_yield_data[index[i]:index[i+1]], xerr=rho_d_2H_WPhi_Phi_width[index[i]:index[i+1]], yerr=rho_d_2H_WPhi_yield_data_statserr[index[i]:index[i+1]], fmt='k.', label='This work')
#     axs[i].set_title(r'$%.1f<\mathrm{E}_{\gamma}<%.1f\ \mathrm{GeV},\  %.1f<-t<%.1f\ \mathrm{GeV}^2$' % (rho_d_2H_WPhi_energy_low[index[i]], rho_d_2H_WPhi_energy_high[index[i]], rho_d_2H_WPhi_minust_low[index[i]], rho_d_2H_WPhi_minust_high[index[i]]))
# axs[0].set_xlim(-180, 180)
# axs[0].set_ylim(0, 150)
# axs[0].set_xticks(np.arange(-120, 180, 60))
# axs[0].set_ylabel(r'$\mathrm{Yield}$')
# fig.suptitle(r"$d(\gamma, \rho d')$ yield vs $\Phi$")
# fig.supxlabel(r'$\Phi\ [\mathrm{deg}]$')
# plt.legend()
# plt.savefig('output/fig_rho_d_2H_WPhi_yield.png', dpi=300)
# plt.close()

# # Plot the efficiency
# fig = plt.figure(figsize=(6*(len(index) - 1), 6))
# gs = fig.add_gridspec(1, len(index) - 1, wspace=0)
# axs = gs.subplots(sharex=True, sharey=True)
# for i in range(len(index) - 1):
#     axs[i].errorbar(rho_d_2H_WPhi_Phi_center[index[i]:index[i+1]], rho_d_2H_WPhi_efficiency[index[i]:index[i+1]], xerr=rho_d_2H_WPhi_Phi_width[index[i]:index[i+1]], yerr=rho_d_2H_WPhi_efficiency_statserr[index[i]:index[i+1]], fmt='k.', label='This work')
#     axs[i].set_title(r'$%.1f<\mathrm{E}_{\gamma}<%.1f\ \mathrm{GeV},\  %.1f<-t<%.1f\ \mathrm{GeV}^2$' % (rho_d_2H_WPhi_energy_low[index[i]], rho_d_2H_WPhi_energy_high[index[i]], rho_d_2H_WPhi_minust_low[index[i]], rho_d_2H_WPhi_minust_high[index[i]]))
# axs[0].set_xlim(-180, 180)
# axs[0].set_ylim(0, 0.1)
# axs[0].set_xticks(np.arange(-120, 180, 60))
# axs[0].set_ylabel(r'$\mathrm{Efficiency}$')
# fig.suptitle(r"$d(\gamma, \rho d')$ efficiency vs $\Phi$")
# fig.supxlabel(r'$\Phi\ [\mathrm{deg}]$')
# plt.legend()
# plt.savefig('output/fig_rho_d_2H_WPhi_efficiency.png', dpi=300)
# plt.close()

# # Plot the results
# fig = plt.figure(figsize=(6*(len(index) - 1), 6))
# gs = fig.add_gridspec(1, len(index) - 1, wspace=0)
# axs = gs.subplots(sharex=True, sharey=True)
# for i in range(len(index) - 1):
#     axs[i].errorbar(rho_d_2H_WPhi_Phi_center[index[i]:index[i+1]], rho_d_2H_WPhi_results[index[i]:index[i+1]], xerr=rho_d_2H_WPhi_Phi_width[index[i]:index[i+1]], yerr=rho_d_2H_WPhi_results_statserr[index[i]:index[i+1]], fmt='k.', label='This work')
#     axs[i].set_title(r'$%.1f<\mathrm{E}_{\gamma}<%.1f\ \mathrm{GeV},\  %.1f<-t<%.1f\ \mathrm{GeV}^2$' % (rho_d_2H_WPhi_energy_low[index[i]], rho_d_2H_WPhi_energy_high[index[i]], rho_d_2H_WPhi_minust_low[index[i]], rho_d_2H_WPhi_minust_high[index[i]]))
#     curve_fit_params, curve_fit_cov = curve_fit(Wphi_func, rho_d_2H_WPhi_Phi_center[index[i]:index[i+1]], rho_d_2H_WPhi_results[index[i]:index[i+1]], p0=[0.0])
#     curve_fit_residuals = rho_d_2H_WPhi_results[index[i]:index[i+1]] - Wphi_func(rho_d_2H_WPhi_Phi_center[index[i]:index[i+1]], curve_fit_params[0])
#     reduced_chi2 = np.sum((curve_fit_residuals/rho_d_2H_WPhi_results_statserr[index[i]:index[i+1]])**2)/(len(rho_d_2H_WPhi_results[index[i]:index[i+1]])-1)
#     axs[i].plot(np.linspace(-180, 180, 360), Wphi_func(np.linspace(-180, 180, 360), curve_fit_params[0]), 'b--', label='Fit')
#     axs[i].text(-175, 1.95, r'$\rho^1_{00}=%.2f\pm%.2f$' % (2*curve_fit_params[0]/3/0.3, np.sqrt(curve_fit_cov[0])), fontsize=10, color='b', ha='left', va='top')
#     axs[i].text(-175, 1.85, r'$\chi^2/\mathrm{dof}=%.2f$' % reduced_chi2, fontsize=10, color='b', ha='left', va='top')
#     axs[i].plot(np.linspace(-180, 180, 360), np.ones(360), 'r--', label='SCHC+NPE')
# axs[0].set_xlim(-180, 180)
# axs[0].set_ylim(0, 2)
# axs[0].set_xticks(np.arange(-120, 180, 60))
# axs[0].set_ylabel(r'$2\pi W(\Phi)$')
# fig.suptitle(r"$d(\gamma, \rho d')$ normalized distribution of $\Phi$")
# fig.supxlabel(r'$\Phi\ [\mathrm{deg}]$')
# plt.legend()
# plt.savefig('output/fig_rho_d_2H_WPhi_results.png', dpi=300)
# plt.close()

# #======================================================================rho_d_2H_Wpsi======================================================================

# # Read the bin edges
# rho_d_2H_Wpsi_energy_low            = np.loadtxt('configs/bins_rho_d_2H_Wpsi.txt')[:,0]
# rho_d_2H_Wpsi_energy_high           = np.loadtxt('configs/bins_rho_d_2H_Wpsi.txt')[:,1]
# rho_d_2H_Wpsi_minust_low            = np.loadtxt('configs/bins_rho_d_2H_Wpsi.txt')[:,2]
# rho_d_2H_Wpsi_minust_high           = np.loadtxt('configs/bins_rho_d_2H_Wpsi.txt')[:,3]
# rho_d_2H_Wpsi_psi_low               = np.loadtxt('configs/bins_rho_d_2H_Wpsi.txt')[:,4]
# rho_d_2H_Wpsi_psi_high              = np.loadtxt('configs/bins_rho_d_2H_Wpsi.txt')[:,5]
# rho_d_2H_Wpsi_psi_center            = (rho_d_2H_Wpsi_psi_high + rho_d_2H_Wpsi_psi_low) / 2
# rho_d_2H_Wpsi_psi_width             = (rho_d_2H_Wpsi_psi_high - rho_d_2H_Wpsi_psi_low) / 2

# # Read the yield numbers
# rho_d_2H_Wpsi_yield_data            = np.loadtxt('output/yield_rho_d_recon_data_2H_exc_Wpsi.txt')[:,6]
# rho_d_2H_Wpsi_yield_sim             = np.loadtxt('output/yield_rho_d_recon_sim_2H_exc_Wpsi.txt')[:,6]
# rho_d_2H_Wpsi_yield_tagged          = np.loadtxt('output/yield_rho_d_thrown_tagged_2H_Wpsi.txt')[:,6]
# rho_d_2H_Wpsi_yield_data_statserr   = np.sqrt(rho_d_2H_Wpsi_yield_data)
# rho_d_2H_Wpsi_yield_sim_statserr    = np.sqrt(rho_d_2H_Wpsi_yield_sim)
# rho_d_2H_Wpsi_yield_tagged_statserr = np.sqrt(rho_d_2H_Wpsi_yield_tagged)

# # Calculate the efficiency
# rho_d_2H_Wpsi_efficiency            = rho_d_2H_Wpsi_yield_sim/rho_d_2H_Wpsi_yield_tagged
# rho_d_2H_Wpsi_efficiency_statserr   = rho_d_2H_Wpsi_efficiency * np.sqrt(1/rho_d_2H_Wpsi_yield_sim + 1/rho_d_2H_Wpsi_yield_tagged)

# # Calculate the results
# rho_d_2H_Wpsi_results               = rho_d_2H_Wpsi_yield_data/rho_d_2H_Wpsi_efficiency  # raw results
# rho_d_2H_Wpsi_results               = normalize_distribution(rho_d_2H_Wpsi_results, rho_d_2H_Wpsi_energy_low, rho_d_2H_Wpsi_minust_low) # normalize to have the sum equal to 1
# rho_d_2H_Wpsi_results               = 2*np.pi*rho_d_2H_Wpsi_results/((rho_d_2H_Wpsi_psi_high - rho_d_2H_Wpsi_psi_low)/180*np.pi)  # normalize to have the integral equal to 2pi
# rho_d_2H_Wpsi_results_statserr      = rho_d_2H_Wpsi_results/np.sqrt(rho_d_2H_Wpsi_yield_data)

# # Find the indices for the different energy and t bins
# index = []
# for i in range(len(rho_d_2H_Wpsi_results)):
#     if (i == 0):
#         index.append(i)
#     elif (i == len(rho_d_2H_Wpsi_results) - 1):
#         index.append(i+1)
#     else:
#         if (rho_d_2H_Wpsi_energy_low[i] != rho_d_2H_Wpsi_energy_low[i-1]) or (rho_d_2H_Wpsi_minust_low[i] != rho_d_2H_Wpsi_minust_low[i-1]):
#             index.append(i)

# # Plot the data yield
# fig = plt.figure(figsize=(6*(len(index) - 1), 6))
# gs = fig.add_gridspec(1, len(index) - 1, wspace=0)
# axs = gs.subplots(sharex=True, sharey=True)
# for i in range(len(index) - 1):
#     axs[i].errorbar(rho_d_2H_Wpsi_psi_center[index[i]:index[i+1]], rho_d_2H_Wpsi_yield_data[index[i]:index[i+1]], xerr=rho_d_2H_Wpsi_psi_width[index[i]:index[i+1]], yerr=rho_d_2H_Wpsi_yield_data_statserr[index[i]:index[i+1]], fmt='k.', label='This work')
#     axs[i].set_title(r'$%.1f<\mathrm{E}_{\gamma}<%.1f\ \mathrm{GeV},\  %.1f<-t<%.1f\ \mathrm{GeV}^2$' % (rho_d_2H_Wpsi_energy_low[index[i]], rho_d_2H_Wpsi_energy_high[index[i]], rho_d_2H_Wpsi_minust_low[index[i]], rho_d_2H_Wpsi_minust_high[index[i]]))
# axs[0].set_xlim(-180, 180)
# axs[0].set_ylim(0, 150)
# axs[0].set_xticks(np.arange(-120, 180, 60))
# axs[0].set_ylabel(r'$\mathrm{Yield}$')
# fig.suptitle(r"$d(\gamma, \rho d')$ yield vs $\psi$")
# fig.supxlabel(r'$\psi\ [\mathrm{deg}]$')
# plt.legend()
# plt.savefig('output/fig_rho_d_2H_Wpsi_yield.png', dpi=300)
# plt.close()

# # Plot the efficiency
# fig = plt.figure(figsize=(6*(len(index) - 1), 6))
# gs = fig.add_gridspec(1, len(index) - 1, wspace=0)
# axs = gs.subplots(sharex=True, sharey=True)
# for i in range(len(index) - 1):
#     axs[i].errorbar(rho_d_2H_Wpsi_psi_center[index[i]:index[i+1]], rho_d_2H_Wpsi_efficiency[index[i]:index[i+1]], xerr=rho_d_2H_Wpsi_psi_width[index[i]:index[i+1]], yerr=rho_d_2H_Wpsi_efficiency_statserr[index[i]:index[i+1]], fmt='k.', label='This work')
#     axs[i].set_title(r'$%.1f<\mathrm{E}_{\gamma}<%.1f\ \mathrm{GeV},\  %.1f<-t<%.1f\ \mathrm{GeV}^2$' % (rho_d_2H_Wpsi_energy_low[index[i]], rho_d_2H_Wpsi_energy_high[index[i]], rho_d_2H_Wpsi_minust_low[index[i]], rho_d_2H_Wpsi_minust_high[index[i]]))
# axs[0].set_xlim(-180, 180)
# axs[0].set_ylim(0, 0.1)
# axs[0].set_xticks(np.arange(-120, 180, 60))
# axs[0].set_ylabel(r'$\mathrm{Efficiency}$')
# fig.suptitle(r"$d(\gamma, \rho d')$ efficiency vs $\psi$")
# fig.supxlabel(r'$\psi\ [\mathrm{deg}]$')
# plt.legend()
# plt.savefig('output/fig_rho_d_2H_Wpsi_efficiency.png', dpi=300)
# plt.close()

# # Plot the results
# fig = plt.figure(figsize=(6*(len(index) - 1), 6))
# gs = fig.add_gridspec(1, len(index) - 1, wspace=0)
# axs = gs.subplots(sharex=True, sharey=True)
# for i in range(len(index) - 1):
#     axs[i].errorbar(rho_d_2H_Wpsi_psi_center[index[i]:index[i+1]], rho_d_2H_Wpsi_results[index[i]:index[i+1]], xerr=rho_d_2H_Wpsi_psi_width[index[i]:index[i+1]], yerr=rho_d_2H_Wpsi_results_statserr[index[i]:index[i+1]], fmt='k.', label='This work')
#     axs[i].set_title(r'$%.1f<\mathrm{E}_{\gamma}<%.1f\ \mathrm{GeV},\  %.1f<-t<%.1f\ \mathrm{GeV}^2$' % (rho_d_2H_Wpsi_energy_low[index[i]], rho_d_2H_Wpsi_energy_high[index[i]], rho_d_2H_Wpsi_minust_low[index[i]], rho_d_2H_Wpsi_minust_high[index[i]]))
#     curve_fit_params, curve_fit_cov = curve_fit(Wphi_func, rho_d_2H_Wpsi_psi_center[index[i]:index[i+1]], rho_d_2H_Wpsi_results[index[i]:index[i+1]], p0=[0.0])
#     curve_fit_residuals = rho_d_2H_Wpsi_results[index[i]:index[i+1]] - Wphi_func(rho_d_2H_Wpsi_psi_center[index[i]:index[i+1]], curve_fit_params[0])
#     reduced_chi2 = np.sum((curve_fit_residuals/rho_d_2H_Wpsi_results_statserr[index[i]:index[i+1]])**2)/(len(rho_d_2H_Wpsi_results[index[i]:index[i+1]])-1)
#     axs[i].plot(np.linspace(-180, 180, 360), Wphi_func(np.linspace(-180, 180, 360), curve_fit_params[0]), 'b--', label='Fit')
#     axs[i].text(-175, 1.95, r'$\rho^1_{1-1}=%.2f\pm%.2f$' % (-curve_fit_params[0]/0.3, np.sqrt(curve_fit_cov[0])), fontsize=10, color='b', ha='left', va='top')
#     axs[i].text(-175, 1.85, r'$\chi^2/\mathrm{dof}=%.2f$' % reduced_chi2, fontsize=10, color='b', ha='left', va='top')
#     axs[i].plot(np.linspace(-180, 180, 360), np.ones(360)+0.3*np.cos(2*np.linspace(-180, 180, 360)/rad_to_deg), 'r--', label='SCHC+NPE')
# axs[0].set_xlim(-180, 180)
# axs[0].set_ylim(0, 2)
# axs[0].set_xticks(np.arange(-120, 180, 60))
# axs[0].set_ylabel(r'$2\pi W(\psi)$')
# fig.suptitle(r"$d(\gamma, \rho d')$ normalized distribution of $\psi$")
# fig.supxlabel(r'$\psi\ [\mathrm{deg}]$')
# plt.legend()
# plt.savefig('output/fig_rho_d_2H_Wpsi_results.png', dpi=300)
# plt.close()

file_pdf.close()