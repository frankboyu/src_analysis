import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages

rad_to_deg = 180/np.pi
file_pdf = PdfPages("/work/halld2/home/boyu/src_analysis/plot/vm_d/output/plots_vm_d_radius.pdf")
###################################################################### SRC-CT #####################################################################################

# Load ds/dt results
phi_d_2H_dsdt_energy_center     = np.loadtxt('output/table_vm_d_dsdt.txt')[:,0]
phi_d_2H_dsdt_minust_center     = np.loadtxt('output/table_vm_d_dsdt.txt')[:,1]
phi_d_2H_dsdt_results           = np.loadtxt('output/table_vm_d_dsdt.txt')[:,2]
phi_d_2H_dsdt_results_statserr  = np.loadtxt('output/table_vm_d_dsdt.txt')[:,3]
phi_d_2H_dsdt_results_p2perr    = np.loadtxt('output/table_vm_d_dsdt.txt')[:,4]
phi_d_2H_dsdt_results_normerr   = np.loadtxt('output/table_vm_d_dsdt.txt')[:,5]
phi_d_2H_dsdt_results_systerr   = np.loadtxt('output/table_vm_d_dsdt.txt')[:,6]

# Identify indices for different energy bins based on gaps in energy_center
index = []
for i in range(len(phi_d_2H_dsdt_results)):
    if (i == 0):
        index.append(i)
    elif (i == len(phi_d_2H_dsdt_results) - 1):
        index.append(i+1)
    else:
        if (phi_d_2H_dsdt_energy_center[i]-phi_d_2H_dsdt_energy_center[i-1] > 1.0):
            index.append(i)

color_list = ['blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan', 'black', 'yellow', 'magenta', 'teal', 'navy', 'maroon', 'lime', 'coral', 'gold', 'silver']

fitting_indices = np.arange(6, 18)
min_points = 6
minust_array = np.zeros((len(fitting_indices)))
radius_array = np.zeros((len(fitting_indices)))
radius_err_array = np.zeros((len(fitting_indices)))

fig = plt.figure(figsize=(8, 6), dpi=300)

phi_d_2H_dsdt_energy_center_68      = phi_d_2H_dsdt_energy_center   [index[0]:index[1]]
phi_d_2H_dsdt_minust_center_68      = phi_d_2H_dsdt_minust_center   [index[0]:index[1]]
phi_d_2H_dsdt_results_68            = phi_d_2H_dsdt_results         [index[0]:index[1]]
phi_d_2H_dsdt_results_statserr_68   = phi_d_2H_dsdt_results_statserr[index[0]:index[1]]
phi_d_2H_dsdt_results_p2perr_68     = phi_d_2H_dsdt_results_p2perr  [index[0]:index[1]]
phi_d_2H_dsdt_results_normerr_68    = phi_d_2H_dsdt_results_normerr [index[0]:index[1]]
phi_d_2H_dsdt_results_systerr_68    = phi_d_2H_dsdt_results_systerr [index[0]:index[1]]


theory_minust = np.loadtxt('theory/temp.txt')[:,0]
theory_dsdt = np.loadtxt('theory/temp.txt')[:,2]

plt.errorbar(phi_d_2H_dsdt_minust_center_68[1:], phi_d_2H_dsdt_results_68[1:], yerr=phi_d_2H_dsdt_results_statserr_68[1:], fmt='b.', label='Data, 5.8-7.8 GeV')

plt.plot(theory_minust, theory_dsdt, 'r-', label='Theory')

plt.xlabel(r'$-t$ (GeV$^2$)', fontsize=16)
plt.ylabel(r'$\frac{d\sigma}{dt}$ (nb/GeV$^2$)', fontsize=16)
plt.yscale('log')
plt.savefig("output/plots_phi_d_theory_intuition.png")

###################################################################### END #####################################################################################

file_pdf.close()