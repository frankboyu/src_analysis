import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages

rad_to_deg = 180/np.pi
file_pdf = PdfPages("/work/halld2/home/boyu/src_analysis/plot/vm_d/output/plots_vm_d_simweight.pdf")
para_list = []
error_list = []

def dsdt_func(minust, a1, b1, a2, b2):
    return a1*np.exp(-b1*minust) + a2*np.exp(-b2*minust)

def lumi(energy_min, energy_max, length):
    lumi_table = np.loadtxt('/work/halld2/home/boyu/src_analysis/flux/output/2H/lumi_summed_2H.txt')
    length_total = 29.5

    integrated_lumi = np.zeros(energy_min.shape, dtype=float)
    for i in range(len(energy_min)):
        for j in range(len(lumi_table)):
            if (lumi_table[j,3] > energy_min[i]) and (lumi_table[j,3] < energy_max[i]):
                integrated_lumi[i] += lumi_table[j][5]

    return integrated_lumi/length_total * length

###################################################################### DATA YIELD #####################################################################################

# Read the bin edges
phi_d_2H_dsdt_energy_center         = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_dsdt_nominal.txt')[:,0]
phi_d_2H_dsdt_energy_width          = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_dsdt_nominal.txt')[:,1]
phi_d_2H_dsdt_energy_low            = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_dsdt_nominal.txt')[:,2]
phi_d_2H_dsdt_energy_high           = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_dsdt_nominal.txt')[:,3]
phi_d_2H_dsdt_minust_center         = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_dsdt_nominal.txt')[:,4]
phi_d_2H_dsdt_minust_width          = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_dsdt_nominal.txt')[:,5]
phi_d_2H_dsdt_minust_low            = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_dsdt_nominal.txt')[:,6]
phi_d_2H_dsdt_minust_high           = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_dsdt_nominal.txt')[:,7]
phi_d_2H_dsdt_yield_data            = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_dsdt_nominal.txt')[:,8]
phi_d_2H_dsdt_yield_data_statserr   = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_dsdt_nominal.txt')[:,9]

# Find the indices for the different energy and t bins
index = []
for i in range(len(phi_d_2H_dsdt_energy_low)):
    if (i == 0):
        index.append(i)
    elif (i == len(phi_d_2H_dsdt_energy_low) - 1):
        index.append(i+1)
    else:
        if (phi_d_2H_dsdt_energy_low[i] != phi_d_2H_dsdt_energy_low[i-1]):
            index.append(i)

###################################################################### ITERATION 0 #####################################################################################

# Simulation yield numbers
phi_d_2H_dsdt_yield_sim_iter0               = np.loadtxt(f'output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_flat_dsdt_simweight_iter0.txt')[:,8]
phi_d_2H_dsdt_yield_sim_statser_iter0       = np.loadtxt(f'output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_flat_dsdt_simweight_iter0.txt')[:,9]
phi_d_2H_dsdt_yield_tagged_iter0            = np.loadtxt(f'output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_flat_dsdt_simweight_iter0.txt')[:,8]
phi_d_2H_dsdt_yield_tagged_statserr_iter0   = np.loadtxt(f'output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_flat_dsdt_simweight_iter0.txt')[:,9]

# Calculate the efficiency and differential cross section
phi_d_2H_dsdt_efficiency_iter0              = phi_d_2H_dsdt_yield_sim_iter0/phi_d_2H_dsdt_yield_tagged_iter0
phi_d_2H_dsdt_efficiency_statserr_iter0     = phi_d_2H_dsdt_efficiency_iter0*np.sqrt((phi_d_2H_dsdt_yield_sim_statser_iter0/phi_d_2H_dsdt_yield_sim_iter0)**2 + (phi_d_2H_dsdt_yield_tagged_statserr_iter0/phi_d_2H_dsdt_yield_tagged_iter0)**2)
phi_d_2H_dsdt_results_iter0                 = phi_d_2H_dsdt_yield_data/phi_d_2H_dsdt_efficiency_iter0/lumi(phi_d_2H_dsdt_energy_low, phi_d_2H_dsdt_energy_high, 28)/(phi_d_2H_dsdt_minust_high-phi_d_2H_dsdt_minust_low)/0.489/1000
phi_d_2H_dsdt_results_statserr_iter0        = phi_d_2H_dsdt_results_iter0*np.sqrt((phi_d_2H_dsdt_yield_data_statserr/phi_d_2H_dsdt_yield_data)**2 + (phi_d_2H_dsdt_efficiency_statserr_iter0/phi_d_2H_dsdt_efficiency_iter0)**2)

# Plot the cross section
fig = plt.figure(figsize=(8, 6), dpi=300)
color_code = ['b', 'k', 'r']
plt.errorbar(phi_d_2H_dsdt_minust_center[index[0]:index[1]],       phi_d_2H_dsdt_results_iter0[index[0]:index[1]],          xerr=phi_d_2H_dsdt_minust_width[index[0]:index[1]],        yerr=phi_d_2H_dsdt_results_statserr_iter0[index[0]:index[1]],            fmt='b.', label='6-8 GeV')
plt.errorbar(phi_d_2H_dsdt_minust_center[index[1]:index[2]],       phi_d_2H_dsdt_results_iter0[index[1]:index[2]],          xerr=phi_d_2H_dsdt_minust_width[index[1]:index[2]],        yerr=phi_d_2H_dsdt_results_statserr_iter0[index[1]:index[2]],            fmt='k.', label='8-9 GeV')
plt.errorbar(phi_d_2H_dsdt_minust_center[index[2]:index[3]],       phi_d_2H_dsdt_results_iter0[index[2]:index[3]],          xerr=phi_d_2H_dsdt_minust_width[index[2]:index[3]],        yerr=phi_d_2H_dsdt_results_statserr_iter0[index[2]:index[3]],            fmt='r.', label='9-11 GeV')

# Fit the cross section with the function
fit_indices = np.where(phi_d_2H_dsdt_minust_center > 0.24)[0]
curve_fit_params, curve_fit_cov = curve_fit(dsdt_func, \
                                            phi_d_2H_dsdt_minust_center[fit_indices], \
                                            phi_d_2H_dsdt_results_iter0[fit_indices], \
                                            sigma=phi_d_2H_dsdt_results_statserr_iter0[fit_indices], \
                                            absolute_sigma=True, p0=[3000, 15, 15, 3])
curve_fit_residuals             = phi_d_2H_dsdt_results_iter0[fit_indices] - dsdt_func(phi_d_2H_dsdt_minust_center[fit_indices], curve_fit_params[0], curve_fit_params[1], curve_fit_params[2], curve_fit_params[3])
reduced_chi2                    = np.sum((curve_fit_residuals/phi_d_2H_dsdt_results_statserr_iter0[fit_indices])**2)/(len(phi_d_2H_dsdt_results_iter0[fit_indices])-4)
plt.plot(np.linspace(0, 2, 100), dsdt_func(np.linspace(0, 2, 100), curve_fit_params[0], curve_fit_params[1], curve_fit_params[2], curve_fit_params[3]), '--', color = 'g', label='Combined fit')
plt.text(0.01, 1, r'$a_1=%.3f\pm%.3f$' % (curve_fit_params[0], np.sqrt(curve_fit_cov[0][0])), fontsize=10, color='b', ha='left', va='top')
plt.text(0.01, 0.7, r'$b_1=%.3f\pm%.3f$' % (curve_fit_params[1], np.sqrt(curve_fit_cov[1][1])), fontsize=10, color='b', ha='left', va='top')
plt.text(0.01, 0.5, r'$a_2=%.3f\pm%.3f$' % (curve_fit_params[2], np.sqrt(curve_fit_cov[2][2])), fontsize=10, color='b', ha='left', va='top')
plt.text(0.01, 0.35, r'$b_2=%.3f\pm%.3f$' % (curve_fit_params[3], np.sqrt(curve_fit_cov[3][3])), fontsize=10, color='b', ha='left', va='top')
plt.text(0.01, 0.25, r'$\chi^2/ndf=%.2f$' % (reduced_chi2), fontsize=10, color='b', ha='left', va='top')

para_list.append(curve_fit_params)
error_list.append(np.sqrt(np.diag(curve_fit_cov)))

# Format the plot
plt.title("Simulation weight iteration 0")
plt.xlabel(r'$-t\ [GeV^2/c]$')
plt.ylabel(r'$d\sigma/dt\ [nb/(GeV^2/c)]$')
plt.xlim(0, 2)
plt.ylim(1e-1, 1e3)
plt.yscale('log')
plt.legend()
file_pdf.savefig()
plt.close()

###################################################################### ITERATION 1 #####################################################################################

# Simulation yield numbers
phi_d_2H_dsdt_yield_sim_iter1               = np.loadtxt(f'output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_flat_dsdt_simweight_iter1.txt')[:,8]
phi_d_2H_dsdt_yield_sim_statser_iter1       = np.loadtxt(f'output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_flat_dsdt_simweight_iter1.txt')[:,9]
phi_d_2H_dsdt_yield_tagged_iter1            = np.loadtxt(f'output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_flat_dsdt_simweight_iter1.txt')[:,8]
phi_d_2H_dsdt_yield_tagged_statserr_iter1   = np.loadtxt(f'output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_flat_dsdt_simweight_iter1.txt')[:,9]

# Calculate the efficiency and differential cross section
phi_d_2H_dsdt_efficiency_iter1              = phi_d_2H_dsdt_yield_sim_iter1/phi_d_2H_dsdt_yield_tagged_iter1
phi_d_2H_dsdt_efficiency_statserr_iter1     = phi_d_2H_dsdt_efficiency_iter1*np.sqrt((phi_d_2H_dsdt_yield_sim_statser_iter1/phi_d_2H_dsdt_yield_sim_iter1)**2 + (phi_d_2H_dsdt_yield_tagged_statserr_iter1/phi_d_2H_dsdt_yield_tagged_iter1)**2)
phi_d_2H_dsdt_results_iter1                 = phi_d_2H_dsdt_yield_data/phi_d_2H_dsdt_efficiency_iter1/lumi(phi_d_2H_dsdt_energy_low, phi_d_2H_dsdt_energy_high, 28)/(phi_d_2H_dsdt_minust_high-phi_d_2H_dsdt_minust_low)/0.489/1000
phi_d_2H_dsdt_results_statserr_iter1        = phi_d_2H_dsdt_results_iter1*np.sqrt((phi_d_2H_dsdt_yield_data_statserr/phi_d_2H_dsdt_yield_data)**2 + (phi_d_2H_dsdt_efficiency_statserr_iter1/phi_d_2H_dsdt_efficiency_iter1)**2)

# Plot the cross section
fig = plt.figure(figsize=(8, 6), dpi=300)
color_code = ['b', 'k', 'r']
plt.errorbar(phi_d_2H_dsdt_minust_center[index[0]:index[1]],       phi_d_2H_dsdt_results_iter1[index[0]:index[1]],          xerr=phi_d_2H_dsdt_minust_width[index[0]:index[1]],        yerr=phi_d_2H_dsdt_results_statserr_iter1[index[0]:index[1]],            fmt='b.', label='6-8 GeV')
plt.errorbar(phi_d_2H_dsdt_minust_center[index[1]:index[2]],       phi_d_2H_dsdt_results_iter1[index[1]:index[2]],          xerr=phi_d_2H_dsdt_minust_width[index[1]:index[2]],        yerr=phi_d_2H_dsdt_results_statserr_iter1[index[1]:index[2]],            fmt='k.', label='8-9 GeV')
plt.errorbar(phi_d_2H_dsdt_minust_center[index[2]:index[3]],       phi_d_2H_dsdt_results_iter1[index[2]:index[3]],          xerr=phi_d_2H_dsdt_minust_width[index[2]:index[3]],        yerr=phi_d_2H_dsdt_results_statserr_iter1[index[2]:index[3]],            fmt='r.', label='9-11 GeV')

# Fit the cross section with the function
fit_indices = np.where(phi_d_2H_dsdt_minust_center > 0.24)[0]
curve_fit_params, curve_fit_cov = curve_fit(dsdt_func, \
                                            phi_d_2H_dsdt_minust_center[fit_indices], \
                                            phi_d_2H_dsdt_results_iter1[fit_indices], \
                                            sigma=phi_d_2H_dsdt_results_statserr_iter1[fit_indices], \
                                            absolute_sigma=True, p0=[3000, 15, 15, 3])
curve_fit_residuals             = phi_d_2H_dsdt_results_iter1[fit_indices] - dsdt_func(phi_d_2H_dsdt_minust_center[fit_indices], curve_fit_params[0], curve_fit_params[1], curve_fit_params[2], curve_fit_params[3])
reduced_chi2                    = np.sum((curve_fit_residuals/phi_d_2H_dsdt_results_statserr_iter1[fit_indices])**2)/(len(phi_d_2H_dsdt_results_iter1[fit_indices])-4)
plt.plot(np.linspace(0, 2, 100), dsdt_func(np.linspace(0, 2, 100), curve_fit_params[0], curve_fit_params[1], curve_fit_params[2], curve_fit_params[3]), '--', color = 'g', label='Combined fit')
plt.text(0.01, 1, r'$a_1=%.3f\pm%.3f$' % (curve_fit_params[0], np.sqrt(curve_fit_cov[0][0])), fontsize=10, color='b', ha='left', va='top')
plt.text(0.01, 0.7, r'$b_1=%.3f\pm%.3f$' % (curve_fit_params[1], np.sqrt(curve_fit_cov[1][1])), fontsize=10, color='b', ha='left', va='top')
plt.text(0.01, 0.5, r'$a_2=%.3f\pm%.3f$' % (curve_fit_params[2], np.sqrt(curve_fit_cov[2][2])), fontsize=10, color='b', ha='left', va='top')
plt.text(0.01, 0.35, r'$b_2=%.3f\pm%.3f$' % (curve_fit_params[3], np.sqrt(curve_fit_cov[3][3])), fontsize=10, color='b', ha='left', va='top')
plt.text(0.01, 0.25, r'$\chi^2/ndf=%.2f$' % (reduced_chi2), fontsize=10, color='b', ha='left', va='top')

para_list.append(curve_fit_params)
error_list.append(np.sqrt(np.diag(curve_fit_cov)))

# Format the plot
plt.title("Simulation weight iteration 1")
plt.xlabel(r'$-t\ [GeV^2/c]$')
plt.ylabel(r'$d\sigma/dt\ [nb/(GeV^2/c)]$')
plt.xlim(0, 2)
plt.ylim(1e-1, 1e3)
plt.yscale('log')
plt.legend()
file_pdf.savefig()
plt.close()

###################################################################### ITERATION 2 #####################################################################################

# Simulation yield numbers
phi_d_2H_dsdt_yield_sim_iter2               = np.loadtxt(f'output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_flat_dsdt_simweight_iter2.txt')[:,8]
phi_d_2H_dsdt_yield_sim_statser_iter2       = np.loadtxt(f'output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_flat_dsdt_simweight_iter2.txt')[:,9]
phi_d_2H_dsdt_yield_tagged_iter2            = np.loadtxt(f'output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_flat_dsdt_simweight_iter2.txt')[:,8]
phi_d_2H_dsdt_yield_tagged_statserr_iter2   = np.loadtxt(f'output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_flat_dsdt_simweight_iter2.txt')[:,9]

# Calculate the efficiency and differential cross section
phi_d_2H_dsdt_efficiency_iter2              = phi_d_2H_dsdt_yield_sim_iter2/phi_d_2H_dsdt_yield_tagged_iter2
phi_d_2H_dsdt_efficiency_statserr_iter2     = phi_d_2H_dsdt_efficiency_iter2*np.sqrt((phi_d_2H_dsdt_yield_sim_statser_iter2/phi_d_2H_dsdt_yield_sim_iter2)**2 + (phi_d_2H_dsdt_yield_tagged_statserr_iter2/phi_d_2H_dsdt_yield_tagged_iter2)**2)
phi_d_2H_dsdt_results_iter2                 = phi_d_2H_dsdt_yield_data/phi_d_2H_dsdt_efficiency_iter2/lumi(phi_d_2H_dsdt_energy_low, phi_d_2H_dsdt_energy_high, 28)/(phi_d_2H_dsdt_minust_high-phi_d_2H_dsdt_minust_low)/0.489/1000
phi_d_2H_dsdt_results_statserr_iter2        = phi_d_2H_dsdt_results_iter2*np.sqrt((phi_d_2H_dsdt_yield_data_statserr/phi_d_2H_dsdt_yield_data)**2 + (phi_d_2H_dsdt_efficiency_statserr_iter2/phi_d_2H_dsdt_efficiency_iter2)**2)

# Plot the cross section
fig = plt.figure(figsize=(8, 6), dpi=300)
color_code = ['b', 'k', 'r']
plt.errorbar(phi_d_2H_dsdt_minust_center[index[0]:index[1]],       phi_d_2H_dsdt_results_iter2[index[0]:index[1]],          xerr=phi_d_2H_dsdt_minust_width[index[0]:index[1]],        yerr=phi_d_2H_dsdt_results_statserr_iter2[index[0]:index[1]],            fmt='b.', label='6-8 GeV')
plt.errorbar(phi_d_2H_dsdt_minust_center[index[1]:index[2]],       phi_d_2H_dsdt_results_iter2[index[1]:index[2]],          xerr=phi_d_2H_dsdt_minust_width[index[1]:index[2]],        yerr=phi_d_2H_dsdt_results_statserr_iter2[index[1]:index[2]],            fmt='k.', label='8-9 GeV')
plt.errorbar(phi_d_2H_dsdt_minust_center[index[2]:index[3]],       phi_d_2H_dsdt_results_iter2[index[2]:index[3]],          xerr=phi_d_2H_dsdt_minust_width[index[2]:index[3]],        yerr=phi_d_2H_dsdt_results_statserr_iter2[index[2]:index[3]],            fmt='r.', label='9-11 GeV')

# Fit the cross section with the function
fit_indices = np.where(phi_d_2H_dsdt_minust_center > 0.24)[0]
curve_fit_params, curve_fit_cov = curve_fit(dsdt_func, \
                                            phi_d_2H_dsdt_minust_center[fit_indices], \
                                            phi_d_2H_dsdt_results_iter2[fit_indices], \
                                            sigma=phi_d_2H_dsdt_results_statserr_iter2[fit_indices], \
                                            absolute_sigma=True, p0=[3000, 15, 15, 3])
curve_fit_residuals             = phi_d_2H_dsdt_results_iter2[fit_indices] - dsdt_func(phi_d_2H_dsdt_minust_center[fit_indices], curve_fit_params[0], curve_fit_params[1], curve_fit_params[2], curve_fit_params[3])
reduced_chi2                    = np.sum((curve_fit_residuals/phi_d_2H_dsdt_results_statserr_iter2[fit_indices])**2)/(len(phi_d_2H_dsdt_results_iter2[fit_indices])-4)
plt.plot(np.linspace(0, 2, 100), dsdt_func(np.linspace(0, 2, 100), curve_fit_params[0], curve_fit_params[1], curve_fit_params[2], curve_fit_params[3]), '--', color = 'g', label='Combined fit')
plt.text(0.01, 1, r'$a_1=%.3f\pm%.3f$' % (curve_fit_params[0], np.sqrt(curve_fit_cov[0][0])), fontsize=10, color='b', ha='left', va='top')
plt.text(0.01, 0.7, r'$b_1=%.3f\pm%.3f$' % (curve_fit_params[1], np.sqrt(curve_fit_cov[1][1])), fontsize=10, color='b', ha='left', va='top')
plt.text(0.01, 0.5, r'$a_2=%.3f\pm%.3f$' % (curve_fit_params[2], np.sqrt(curve_fit_cov[2][2])), fontsize=10, color='b', ha='left', va='top')
plt.text(0.01, 0.35, r'$b_2=%.3f\pm%.3f$' % (curve_fit_params[3], np.sqrt(curve_fit_cov[3][3])), fontsize=10, color='b', ha='left', va='top')
plt.text(0.01, 0.25, r'$\chi^2/ndf=%.2f$' % (reduced_chi2), fontsize=10, color='b', ha='left', va='top')

para_list.append(curve_fit_params)
error_list.append(np.sqrt(np.diag(curve_fit_cov)))

# Format the plot
plt.title("Simulation weight iteration 2")
plt.xlabel(r'$-t\ [GeV^2/c]$')
plt.ylabel(r'$d\sigma/dt\ [nb/(GeV^2/c)]$')
plt.xlim(0, 2)
plt.ylim(1e-1, 1e3)
plt.yscale('log')
plt.legend()
file_pdf.savefig()
plt.close()

###################################################################### ITERATION 3 #####################################################################################

# Simulation yield numbers
phi_d_2H_dsdt_yield_sim_iter3               = np.loadtxt(f'output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_flat_dsdt_simweight_iter3.txt')[:,8]
phi_d_2H_dsdt_yield_sim_statser_iter3       = np.loadtxt(f'output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_flat_dsdt_simweight_iter3.txt')[:,9]
phi_d_2H_dsdt_yield_tagged_iter3            = np.loadtxt(f'output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_flat_dsdt_simweight_iter3.txt')[:,8]
phi_d_2H_dsdt_yield_tagged_statserr_iter3   = np.loadtxt(f'output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_flat_dsdt_simweight_iter3.txt')[:,9]

# Calculate the efficiency and differential cross section
phi_d_2H_dsdt_efficiency_iter3              = phi_d_2H_dsdt_yield_sim_iter3/phi_d_2H_dsdt_yield_tagged_iter3
phi_d_2H_dsdt_efficiency_statserr_iter3     = phi_d_2H_dsdt_efficiency_iter3*np.sqrt((phi_d_2H_dsdt_yield_sim_statser_iter3/phi_d_2H_dsdt_yield_sim_iter3)**2 + (phi_d_2H_dsdt_yield_tagged_statserr_iter3/phi_d_2H_dsdt_yield_tagged_iter3)**2)
phi_d_2H_dsdt_results_iter3                 = phi_d_2H_dsdt_yield_data/phi_d_2H_dsdt_efficiency_iter3/lumi(phi_d_2H_dsdt_energy_low, phi_d_2H_dsdt_energy_high, 28)/(phi_d_2H_dsdt_minust_high-phi_d_2H_dsdt_minust_low)/0.489/1000
phi_d_2H_dsdt_results_statserr_iter3        = phi_d_2H_dsdt_results_iter3*np.sqrt((phi_d_2H_dsdt_yield_data_statserr/phi_d_2H_dsdt_yield_data)**2 + (phi_d_2H_dsdt_efficiency_statserr_iter3/phi_d_2H_dsdt_efficiency_iter3)**2)

# Plot the cross section
fig = plt.figure(figsize=(8, 6), dpi=300)
color_code = ['b', 'k', 'r']
plt.errorbar(phi_d_2H_dsdt_minust_center[index[0]:index[1]],       phi_d_2H_dsdt_results_iter3[index[0]:index[1]],          xerr=phi_d_2H_dsdt_minust_width[index[0]:index[1]],        yerr=phi_d_2H_dsdt_results_statserr_iter3[index[0]:index[1]],            fmt='b.', label='6-8 GeV')
plt.errorbar(phi_d_2H_dsdt_minust_center[index[1]:index[2]],       phi_d_2H_dsdt_results_iter3[index[1]:index[2]],          xerr=phi_d_2H_dsdt_minust_width[index[1]:index[2]],        yerr=phi_d_2H_dsdt_results_statserr_iter3[index[1]:index[2]],            fmt='k.', label='8-9 GeV')
plt.errorbar(phi_d_2H_dsdt_minust_center[index[2]:index[3]],       phi_d_2H_dsdt_results_iter3[index[2]:index[3]],          xerr=phi_d_2H_dsdt_minust_width[index[2]:index[3]],        yerr=phi_d_2H_dsdt_results_statserr_iter3[index[2]:index[3]],            fmt='r.', label='9-11 GeV')

# Fit the cross section with the function
fit_indices = np.where(phi_d_2H_dsdt_minust_center > 0.24)[0]
curve_fit_params, curve_fit_cov = curve_fit(dsdt_func, \
                                            phi_d_2H_dsdt_minust_center[fit_indices], \
                                            phi_d_2H_dsdt_results_iter3[fit_indices], \
                                            sigma=phi_d_2H_dsdt_results_statserr_iter3[fit_indices], \
                                            absolute_sigma=True, p0=[3000, 15, 15, 3])
curve_fit_residuals             = phi_d_2H_dsdt_results_iter3[fit_indices] - dsdt_func(phi_d_2H_dsdt_minust_center[fit_indices], curve_fit_params[0], curve_fit_params[1], curve_fit_params[2], curve_fit_params[3])
reduced_chi2                    = np.sum((curve_fit_residuals/phi_d_2H_dsdt_results_statserr_iter3[fit_indices])**2)/(len(phi_d_2H_dsdt_results_iter3[fit_indices])-4)
plt.plot(np.linspace(0, 2, 100), dsdt_func(np.linspace(0, 2, 100), curve_fit_params[0], curve_fit_params[1], curve_fit_params[2], curve_fit_params[3]), '--', color = 'g', label='Combined fit')
plt.text(0.01, 1, r'$a_1=%.3f\pm%.3f$' % (curve_fit_params[0], np.sqrt(curve_fit_cov[0][0])), fontsize=10, color='b', ha='left', va='top')
plt.text(0.01, 0.7, r'$b_1=%.3f\pm%.3f$' % (curve_fit_params[1], np.sqrt(curve_fit_cov[1][1])), fontsize=10, color='b', ha='left', va='top')
plt.text(0.01, 0.5, r'$a_2=%.3f\pm%.3f$' % (curve_fit_params[2], np.sqrt(curve_fit_cov[2][2])), fontsize=10, color='b', ha='left', va='top')
plt.text(0.01, 0.35, r'$b_2=%.3f\pm%.3f$' % (curve_fit_params[3], np.sqrt(curve_fit_cov[3][3])), fontsize=10, color='b', ha='left', va='top')
plt.text(0.01, 0.25, r'$\chi^2/ndf=%.2f$' % (reduced_chi2), fontsize=10, color='b', ha='left', va='top')

para_list.append(curve_fit_params)
error_list.append(np.sqrt(np.diag(curve_fit_cov)))

# Format the plot
plt.title("Simulation weight iteration 3")
plt.xlabel(r'$-t\ [GeV^2/c]$')
plt.ylabel(r'$d\sigma/dt\ [nb/(GeV^2/c)]$')
plt.xlim(0, 2)
plt.ylim(1e-1, 1e3)
plt.yscale('log')
plt.legend()
file_pdf.savefig()
plt.close()

###################################################################### SUMMARY #####################################################################################


# a1_list = np.array([3632.01, 6127.78, 6417.46, 6442.44, 6444.93, 6444.97, 6445.03])
# b1_list = np.array([15.82, 17.10, 17.21, 17.22, 17.22, 17.22, 17.22])
# a2_list = np.array([12.19, 12.43, 12.42, 12.42, 12.42, 12.42, 12.42])
# b2_list = np.array([2.49, 2.50, 2.49, 2.49, 2.49, 2.49, 2.49])

# # plt.plot(np.arange(len(a1_list)), (a1_list-a1_list[-1])/a1_list[-1], 'bo-', label='a1')
# # plt.plot(np.arange(len(b1_list)), (b1_list-b1_list[-1])/b1_list[-1], 'ro-', label='b1')
# # plt.plot(np.arange(len(a2_list)), (a2_list-a2_list[-1])/a2_list[-1], 'go-', label='a2')
# # plt.plot(np.arange(len(b2_list)), (b2_list-b2_list[-1])/b2_list[-1], 'mo-', label='b2')
# plt.plot(np.arange(1, len(a1_list)), (a1_list[0:-1]-a1_list[1:])/a1_list[0:-1], 'bo-', label='a1')
# plt.plot(np.arange(1, len(b1_list)), (b1_list[0:-1]-b1_list[1:])/b1_list[0:-1], 'ro-', label='b1')
# plt.plot(np.arange(1, len(a2_list)), (a2_list[0:-1]-a2_list[1:])/a2_list[0:-1], 'go-', label='a2')
# plt.plot(np.arange(1, len(b2_list)), (b2_list[0:-1]-b2_list[1:])/b2_list[0:-1], 'mo-', label='b2')
# plt.axhline(y=0, color='k', linestyle='--')
# plt.title(r"Fit parameter variation with different simulation weights")
# plt.xlabel('Iteration')
# plt.ylabel('Relative variation')
# plt.legend()
# file_pdf.savefig()
# plt.close()

file_pdf.close()