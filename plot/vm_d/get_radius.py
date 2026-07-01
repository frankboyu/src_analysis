import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages

rad_to_deg = 180/np.pi
file_pdf = PdfPages("/work/halld2/home/boyu/src_analysis/plot/vm_d/output/plots_vm_d_radius.pdf")

def monopole_form(minust, norm, masssq):
    form_factor = 1 / (1 + minust/masssq)**1
    return (norm * form_factor)**2

def dipole_form(minust, norm, masssq):
    form_factor = 1 / (1 + minust/masssq)**2
    return (norm * form_factor)**2

def tripole_form(minust, norm, masssq):
    form_factor = 1 / (1 + minust/masssq)**3
    return (norm * form_factor)**2

def quadrupole_form(minust, norm, masssq):
    form_factor = 1 / (1 + minust/masssq)**4
    return (norm * form_factor)**2

def gaussian_form(minust, norm, masssq):
    form_factor = np.exp(-minust/4/masssq)
    return (norm * form_factor)**2

def empirical_form(minust, norm, masssq):
    form_factor = norm * np.exp(-minust/masssq)
    return (norm * form_factor)**2

def uniform_form(minust, norm, radius):
    x = np.sqrt(minust)*radius
    form_factor = np.sin(x) / x**3 - np.cos(x) / x**2
    return (norm * form_factor)**2

def rational_1_1(minust, p0, p1, p2):
    form_factor = (1 + p1*minust) / (1 + p2*minust)
    return (p0 * form_factor)**2

def polynomial_1(minust, p0, p1):
    form_factor = 1 / (1 + p1*minust)
    return (p0 * form_factor)**2

# def continued_fraction_1(minust, p0, p1, p2):
#     p_list = np.array([p0, p1, p2])
#     form_factor = 1
#     for i in range(1):
#         form_factor = 1 + 
#     form_factor = 1 / (1 + p1*minust + p2*minust**2)
#     return (p0 * form_factor)**2

###################################################################### CLAS #####################################################################################
# CLAS results
phi_d_2H_dsdt_clas_minust_low           = np.array([0.350, 0.375, 0.400, 0.425, 0.450, 0.500, 0.550, 0.600, 0.700, 0.800, 1.000, 1.200, 1.400])
phi_d_2H_dsdt_clas_minust_high          = np.array([0.375, 0.400, 0.425, 0.450, 0.500, 0.550, 0.600, 0.700, 0.800, 1.000, 1.200, 1.400, 2.000])
phi_d_2H_dsdt_clas_minust_middle        = (phi_d_2H_dsdt_clas_minust_high + phi_d_2H_dsdt_clas_minust_low) / 2
phi_d_2H_dsdt_clas_minust_size          = (phi_d_2H_dsdt_clas_minust_high - phi_d_2H_dsdt_clas_minust_low) / 2
phi_d_2H_dsdt_clas_minust_center        = np.array([0.360,  0.385,  0.410,  0.435,  0.474,  0.524,  0.574,  0.646,  0.746,  0.888,  1.091,  1.292,  1.637])

phi_d_2H_dsdt_clas_results_16           = np.array([10.21,  8.85,   7.32,   6.16,   4.73,   3.52,   2.66,   2.17,   1.40,   0.94,   0.57,   0.28,   0.19])
phi_d_2H_dsdt_clas_results_16_statserr  = np.array([0.82,   0.75,   0.59,   0.55,   0.34,   0.28,   0.24,   0.15,   0.12,   0.07,   0.06,   0.05,   0.02])
phi_d_2H_dsdt_clas_results_16_systerr   = np.array([1.70,   1.11,   0.94,   0.81,   0.60,   0.51,   0.38,   0.26,   0.16,   0.11,   0.07,   0.04,   0.03])
phi_d_2H_dsdt_clas_results_16_totalerr  = np.sqrt(phi_d_2H_dsdt_clas_results_16_statserr**2 + phi_d_2H_dsdt_clas_results_16_systerr**2)

phi_d_2H_dsdt_clas_results_26           = np.array([8.63,   6.80,   4.57,   5.76,   3.99,   3.59,   2.11,   1.83,   1.32,   0.96,   0.57,   0.36,   0.15])
phi_d_2H_dsdt_clas_results_26_statserr  = np.array([0.80,   0.69,   0.53,   0.56,   0.33,   0.29,   0.22,   0.14,   0.12,   0.07,   0.05,   0.04,   0.02])
phi_d_2H_dsdt_clas_results_26_systerr   = np.array([1.04,   1.07,   0.74,   0.65,   0.55,   0.55,   0.28,   0.24,   0.20,   0.11,   0.06,   0.05,   0.02])
phi_d_2H_dsdt_clas_results_26_totalerr  = np.sqrt(phi_d_2H_dsdt_clas_results_26_statserr**2 + phi_d_2H_dsdt_clas_results_26_systerr**2)

# Plot the cross section
fig = plt.figure(figsize=(8, 6), dpi=300)

plt.errorbar(phi_d_2H_dsdt_clas_minust_middle, phi_d_2H_dsdt_clas_results_16, yerr=phi_d_2H_dsdt_clas_results_16_statserr, fmt='o', color='blue', label='CLAS 16')

fit_indices = np.where((phi_d_2H_dsdt_clas_minust_middle < 0.60) & (phi_d_2H_dsdt_clas_minust_middle > 0.1))[0]

curve_fit_params, curve_fit_cov = curve_fit(monopole_form, \
                                            phi_d_2H_dsdt_clas_minust_middle[fit_indices], \
                                            phi_d_2H_dsdt_clas_results_16[fit_indices], \
                                            sigma=phi_d_2H_dsdt_clas_results_16_statserr[fit_indices],
                                            absolute_sigma=True, p0=[100, 0.05], maxfev=5000)
curve_fit_residuals = phi_d_2H_dsdt_clas_results_16[fit_indices] - monopole_form(phi_d_2H_dsdt_clas_minust_middle[fit_indices], curve_fit_params[0], curve_fit_params[1])
reduced_chi2 = np.sum(curve_fit_residuals**2) / (len(fit_indices) - len(curve_fit_params))
plt.plot(np.linspace(0.001, 0.6, 100), monopole_form(np.linspace(0.001, 0.6, 100), curve_fit_params[0], curve_fit_params[1]), '-', color='blue', label=f'Monopole fit, χ²/ndf = {reduced_chi2:.2f}, r={0.197*np.sqrt(6/curve_fit_params[1]):.2f} fm')

curve_fit_params, curve_fit_cov = curve_fit(dipole_form, \
                                            phi_d_2H_dsdt_clas_minust_middle[fit_indices], \
                                            phi_d_2H_dsdt_clas_results_16[fit_indices], \
                                            sigma=phi_d_2H_dsdt_clas_results_16_statserr[fit_indices],
                                            absolute_sigma=True, p0=[100, 0.05], maxfev=1000)
curve_fit_residuals = phi_d_2H_dsdt_clas_results_16[fit_indices] - dipole_form(phi_d_2H_dsdt_clas_minust_middle[fit_indices], curve_fit_params[0], curve_fit_params[1])
reduced_chi2 = np.sum(curve_fit_residuals**2) / (len(fit_indices) - len(curve_fit_params))
plt.plot(np.linspace(0.001, 0.6, 100), dipole_form(np.linspace(0.001, 0.6, 100), curve_fit_params[0], curve_fit_params[1]), '--', color='black', label=f'Dipole fit, χ²/ndf = {reduced_chi2:.2f}, r={0.1973*np.sqrt(12/curve_fit_params[1]):.2f} fm')

curve_fit_params, curve_fit_cov = curve_fit(tripole_form, \
                                            phi_d_2H_dsdt_clas_minust_middle[fit_indices], \
                                            phi_d_2H_dsdt_clas_results_16[fit_indices], \
                                            sigma=phi_d_2H_dsdt_clas_results_16_statserr[fit_indices],
                                            absolute_sigma=True, p0=[100, 0.05], maxfev=1000)
curve_fit_residuals = phi_d_2H_dsdt_clas_results_16[fit_indices] - tripole_form(phi_d_2H_dsdt_clas_minust_middle[fit_indices], curve_fit_params[0], curve_fit_params[1])
reduced_chi2 = np.sum(curve_fit_residuals**2) / (len(fit_indices) - len(curve_fit_params))
plt.plot(np.linspace(0.001, 0.6, 100), tripole_form(np.linspace(0.001, 0.6, 100), curve_fit_params[0], curve_fit_params[1]), '-.', color='red', label=f'Tripole fit, χ²/ndf = {reduced_chi2:.2f}, r={0.1973*np.sqrt(18/curve_fit_params[1]):.2f} fm')

curve_fit_params, curve_fit_cov = curve_fit(quadrupole_form, \
                                            phi_d_2H_dsdt_clas_minust_middle[fit_indices], \
                                            phi_d_2H_dsdt_clas_results_16[fit_indices], \
                                            sigma=phi_d_2H_dsdt_clas_results_16_statserr[fit_indices],
                                            absolute_sigma=True, p0=[100, 0.05], maxfev=1000)
curve_fit_residuals = phi_d_2H_dsdt_clas_results_16[fit_indices] - quadrupole_form(phi_d_2H_dsdt_clas_minust_middle[fit_indices], curve_fit_params[0], curve_fit_params[1])
reduced_chi2 = np.sum(curve_fit_residuals**2) / (len(fit_indices) - len(curve_fit_params))
plt.plot(np.linspace(0.001, 0.6, 100), quadrupole_form(np.linspace(0.001, 0.6, 100), curve_fit_params[0], curve_fit_params[1]), ':', color='green', label=f'Quadrupole fit, χ²/ndf = {reduced_chi2:.2f}, r={0.1973*np.sqrt(24/curve_fit_params[1]):.2f} fm')

curve_fit_params, curve_fit_cov = curve_fit(gaussian_form, \
                                            phi_d_2H_dsdt_clas_minust_middle[fit_indices], \
                                            phi_d_2H_dsdt_clas_results_16[fit_indices], \
                                            sigma=phi_d_2H_dsdt_clas_results_16_statserr[fit_indices],
                                            absolute_sigma=True, p0=[100, 0.05], maxfev=1000)
curve_fit_residuals = phi_d_2H_dsdt_clas_results_16[fit_indices] - gaussian_form(phi_d_2H_dsdt_clas_minust_middle[fit_indices], curve_fit_params[0], curve_fit_params[1])
reduced_chi2 = np.sum(curve_fit_residuals**2) / (len(fit_indices) - len(curve_fit_params))
plt.plot(np.linspace(0.001, 0.6, 100), gaussian_form(np.linspace(0.001, 0.6, 100), curve_fit_params[0], curve_fit_params[1]), '-.', color='purple', label=f'Gaussian fit, χ²/ndf = {reduced_chi2:.2f}, r={0.1973*np.sqrt(3/2/curve_fit_params[1]):.2f} fm')

curve_fit_params, curve_fit_cov = curve_fit(uniform_form, \
                                            phi_d_2H_dsdt_clas_minust_middle[fit_indices], \
                                            phi_d_2H_dsdt_clas_results_16[fit_indices], \
                                            sigma=phi_d_2H_dsdt_clas_results_16_statserr[fit_indices],
                                            absolute_sigma=True, p0=[100, 0.05], maxfev=1000)
curve_fit_residuals = phi_d_2H_dsdt_clas_results_16[fit_indices] - uniform_form(phi_d_2H_dsdt_clas_minust_middle[fit_indices], curve_fit_params[0], curve_fit_params[1])
reduced_chi2 = np.sum(curve_fit_residuals**2) / (len(fit_indices) - len(curve_fit_params))
plt.plot(np.linspace(0.001, 0.6, 100), uniform_form(np.linspace(0.001, 0.6, 100), curve_fit_params[0], curve_fit_params[1]), '-.', color='orange', label=f'Uniform fit, χ²/ndf = {reduced_chi2:.2f}, r={0.1973*np.sqrt(3*curve_fit_params[1]*curve_fit_params[1]/5):.2f} fm')

# curve_fit_params, curve_fit_cov = curve_fit(rational_1_1, \
#                                             phi_d_2H_dsdt_clas_minust_middle[fit_indices], \
#                                             phi_d_2H_dsdt_clas_results_16[fit_indices], \
#                                             sigma=phi_d_2H_dsdt_clas_results_16_statserr[fit_indices],
#                                             absolute_sigma=True, p0=[100, 0.05, 0.05], maxfev=10000)
# curve_fit_residuals = phi_d_2H_dsdt_clas_results_16[fit_indices] - rational_1_1(phi_d_2H_dsdt_clas_minust_middle[fit_indices], curve_fit_params[0], curve_fit_params[1], curve_fit_params[2])
# reduced_chi2 = np.sum(curve_fit_residuals**2) / (len(fit_indices) - len(curve_fit_params))
# plt.plot(np.linspace(0.001, 0.6, 100), rational_1_1(np.linspace(0.001, 0.6, 100), curve_fit_params[0], curve_fit_params[1], curve_fit_params[2]), '-', color='red', label=f'Rational fit (χ²/ndf = {reduced_chi2:.2f})')


# Format the plot
plt.title("CLAS 1.6-2.6 GeV")
plt.xlabel(r'$-t\ [GeV^2/c]$')
plt.ylabel(r'$d\sigma/dt\ [nb/(GeV^2/c)]$')
plt.xlim(0.0, 0.6)
plt.ylim(2e0, 2e2)
plt.yscale('log')
plt.legend()
file_pdf.savefig()
plt.close()

###################################################################### LEPS #####################################################################################

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

###################################################################### 5.8-7.8 GeV #####################################################################################

phi_d_2H_dsdt_energy_center_68      = phi_d_2H_dsdt_energy_center   [index[0]:index[1]]
phi_d_2H_dsdt_minust_center_68      = phi_d_2H_dsdt_minust_center   [index[0]:index[1]]
phi_d_2H_dsdt_results_68            = phi_d_2H_dsdt_results         [index[0]:index[1]]
phi_d_2H_dsdt_results_statserr_68   = phi_d_2H_dsdt_results_statserr[index[0]:index[1]]
phi_d_2H_dsdt_results_p2perr_68     = phi_d_2H_dsdt_results_p2perr  [index[0]:index[1]]
phi_d_2H_dsdt_results_normerr_68    = phi_d_2H_dsdt_results_normerr [index[0]:index[1]]
phi_d_2H_dsdt_results_systerr_68    = phi_d_2H_dsdt_results_systerr [index[0]:index[1]]

# MONOPOLE FORM
fig = plt.figure(figsize=(16, 6), dpi=300)
gs = fig.add_gridspec(1, 2)
axs = gs.subplots()
plt.suptitle("Monopole Form, 5.8-7.8 GeV")
plt.axes(axs[0])
plt.errorbar(phi_d_2H_dsdt_minust_center_68[1:], phi_d_2H_dsdt_results_68[1:], yerr=phi_d_2H_dsdt_results_statserr_68[1:], fmt='b.', label='Data, 5.8-7.8 GeV')

for i in range(6, 16):
    curve_fit_params, curve_fit_cov = curve_fit(monopole_form, \
                                            phi_d_2H_dsdt_minust_center_68[1:i], \
                                            phi_d_2H_dsdt_results_68[1:i], \
                                            sigma=phi_d_2H_dsdt_results_statserr_68[1:i], \
                                            absolute_sigma=True, p0=[300, 0.05], maxfev=50000)
    curve_fit_residuals             = phi_d_2H_dsdt_results_68[1:i] - monopole_form(phi_d_2H_dsdt_minust_center_68[1:i], curve_fit_params[0], curve_fit_params[1])
    reduced_chi2                    = np.sum((curve_fit_residuals/phi_d_2H_dsdt_results_statserr_68[1:i])**2)/(len(phi_d_2H_dsdt_results_68[1:i])-2)
    plt.axes(axs[0])
    plt.plot(np.linspace(0.01, phi_d_2H_dsdt_minust_center_68[i], 1000), monopole_form(np.linspace(0.01, phi_d_2H_dsdt_minust_center_68[i], 1000), curve_fit_params[0], curve_fit_params[1]), '-',  color = color_list[i-3])
    plt.plot(np.linspace(phi_d_2H_dsdt_minust_center_68[i], 0.6, 1000),  monopole_form(np.linspace(phi_d_2H_dsdt_minust_center_68[i], 0.6, 1000), curve_fit_params[0],  curve_fit_params[1]), ':', color = color_list[i-3])

    plt.axes(axs[1])
    # plt.errorbar(phi_d_2H_dsdt_minust_center_68[i], 0.197*np.sqrt(6/curve_fit_params[1]), yerr=0.197*np.sqrt(6/curve_fit_cov[1][1]), fmt='o', color = color_list[i-3])

plt.axes(axs[0])
plt.yscale('log')
plt.xlabel(r'$-t\ [GeV^2/c]$')
plt.ylabel(r'$d\sigma/dt\ [nb/(GeV^2/c)]$')
plt.xlim(0, 0.6)
plt.ylim(1e0, 1e3)
plt.legend()

plt.axes(axs[1])
plt.xlabel(r'$-t\ [GeV^2/c]$')
plt.ylabel(r'$r\ [fm]$')
plt.xlim(0, 0.6)
plt.ylim(0, 2.5)

file_pdf.savefig()
plt.close()


# DIPOLE FORM
fig = plt.figure(figsize=(16, 6), dpi=300)
gs = fig.add_gridspec(1, 2)
axs = gs.subplots()
plt.suptitle("Dipole Form, 5.8-7.8 GeV")
plt.axes(axs[0])
plt.errorbar(phi_d_2H_dsdt_minust_center_68[1:], phi_d_2H_dsdt_results_68[1:], yerr=phi_d_2H_dsdt_results_statserr_68[1:], fmt='b.', label='Data, 5.8-7.8 GeV')

for i in range(6, 16):
    curve_fit_params, curve_fit_cov = curve_fit(dipole_form, \
                                            phi_d_2H_dsdt_minust_center_68[1:i], \
                                            phi_d_2H_dsdt_results_68[1:i], \
                                            sigma=phi_d_2H_dsdt_results_statserr_68[1:i], \
                                            absolute_sigma=True, p0=[300, 0.05], maxfev=50000)
    curve_fit_residuals             = phi_d_2H_dsdt_results_68[1:i] - dipole_form(phi_d_2H_dsdt_minust_center_68[1:i], curve_fit_params[0], curve_fit_params[1])
    reduced_chi2                    = np.sum((curve_fit_residuals/phi_d_2H_dsdt_results_statserr_68[1:i])**2)/(len(phi_d_2H_dsdt_results_68[1:i])-2)
    plt.axes(axs[0])
    plt.plot(np.linspace(0.01, phi_d_2H_dsdt_minust_center_68[i], 1000), dipole_form(np.linspace(0.01, phi_d_2H_dsdt_minust_center_68[i], 1000), curve_fit_params[0], curve_fit_params[1]), '-',  color = color_list[i-3])
    plt.plot(np.linspace(phi_d_2H_dsdt_minust_center_68[i], 0.6, 1000),  dipole_form(np.linspace(phi_d_2H_dsdt_minust_center_68[i], 0.6, 1000), curve_fit_params[0],  curve_fit_params[1]), ':', color = color_list[i-3])

    plt.axes(axs[1])
    # plt.errorbar(phi_d_2H_dsdt_minust_center_68[i], 0.197*np.sqrt(6/curve_fit_params[1]), yerr=0.197*np.sqrt(6/curve_fit_cov[1][1]), fmt='o', color = color_list[i-3])

plt.axes(axs[0])
plt.yscale('log')
plt.xlabel(r'$-t\ [GeV^2/c]$')
plt.ylabel(r'$d\sigma/dt\ [nb/(GeV^2/c)]$')
plt.xlim(0, 0.6)
plt.ylim(1e0, 1e3)
plt.legend()

plt.axes(axs[1])
plt.xlabel(r'$-t\ [GeV^2/c]$')
plt.ylabel(r'$r\ [fm]$')
plt.xlim(0, 0.6)
plt.ylim(0, 2.5)

file_pdf.savefig()
plt.close()

# TRIPOLE FORM
fig = plt.figure(figsize=(16, 6), dpi=300)
gs = fig.add_gridspec(1, 2)
axs = gs.subplots()
plt.suptitle("Tripole Form, 5.8-7.8 GeV")
plt.axes(axs[0])
plt.errorbar(phi_d_2H_dsdt_minust_center_68[1:], phi_d_2H_dsdt_results_68[1:], yerr=phi_d_2H_dsdt_results_statserr_68[1:], fmt='b.', label='Data, 5.8-7.8 GeV')

for i in range(6, 16):
    curve_fit_params, curve_fit_cov = curve_fit(tripole_form, \
                                            phi_d_2H_dsdt_minust_center_68[1:i], \
                                            phi_d_2H_dsdt_results_68[1:i], \
                                            sigma=phi_d_2H_dsdt_results_statserr_68[1:i], \
                                            absolute_sigma=True, p0=[300, 0.05], maxfev=50000)
    curve_fit_residuals             = phi_d_2H_dsdt_results_68[1:i] - tripole_form(phi_d_2H_dsdt_minust_center_68[1:i], curve_fit_params[0], curve_fit_params[1])
    reduced_chi2                    = np.sum((curve_fit_residuals/phi_d_2H_dsdt_results_statserr_68[1:i])**2)/(len(phi_d_2H_dsdt_results_68[1:i])-2)
    plt.axes(axs[0])
    plt.plot(np.linspace(0.01, phi_d_2H_dsdt_minust_center_68[i], 1000), tripole_form(np.linspace(0.01, phi_d_2H_dsdt_minust_center_68[i], 1000), curve_fit_params[0], curve_fit_params[1]), '-',  color = color_list[i-3])
    plt.plot(np.linspace(phi_d_2H_dsdt_minust_center_68[i], 0.6, 1000),  tripole_form(np.linspace(phi_d_2H_dsdt_minust_center_68[i], 0.6, 1000), curve_fit_params[0],  curve_fit_params[1]), ':', color = color_list[i-3])

    plt.axes(axs[1])
    this_radius = 0.1973*np.sqrt(18/curve_fit_params[1])
    this_radius_err = 0.1973*np.sqrt(18)*0.5*np.sqrt(curve_fit_cov[1][1])/curve_fit_params[1]**(3/2)
    plt.errorbar(phi_d_2H_dsdt_minust_center_68[i], this_radius, yerr=this_radius_err, fmt='o', color = color_list[i-3])

plt.axes(axs[0])
plt.yscale('log')
plt.xlabel(r'$-t\ [GeV^2/c]$')
plt.ylabel(r'$d\sigma/dt\ [nb/(GeV^2/c)]$')
plt.xlim(0, 0.6)
plt.ylim(1e0, 1e3)
plt.legend()

plt.axes(axs[1])
plt.xlabel(r'$-t\ [GeV^2/c]$')
plt.ylabel(r'$r\ [fm]$')
plt.xlim(0, 0.6)
plt.ylim(0, 2.5)

file_pdf.savefig()
plt.close()

# QUADRUPOLE FORM
fig = plt.figure(figsize=(16, 6), dpi=300)
gs = fig.add_gridspec(1, 2)
axs = gs.subplots()
plt.suptitle("Quadrupole Form, 5.8-7.8 GeV")
plt.axes(axs[0])
plt.errorbar(phi_d_2H_dsdt_minust_center_68[1:], phi_d_2H_dsdt_results_68[1:], yerr=phi_d_2H_dsdt_results_statserr_68[1:], fmt='b.', label='Data, 5.8-7.8 GeV')

for i in range(6, 16):
    curve_fit_params, curve_fit_cov = curve_fit(quadrupole_form, \
                                            phi_d_2H_dsdt_minust_center_68[1:i], \
                                            phi_d_2H_dsdt_results_68[1:i], \
                                            sigma=phi_d_2H_dsdt_results_statserr_68[1:i], \
                                            absolute_sigma=True, p0=[300, 0.05], maxfev=50000)
    curve_fit_residuals             = phi_d_2H_dsdt_results_68[1:i] - quadrupole_form(phi_d_2H_dsdt_minust_center_68[1:i], curve_fit_params[0], curve_fit_params[1])
    reduced_chi2                    = np.sum((curve_fit_residuals/phi_d_2H_dsdt_results_statserr_68[1:i])**2)/(len(phi_d_2H_dsdt_results_68[1:i])-2)
    plt.axes(axs[0])
    plt.plot(np.linspace(0.01, phi_d_2H_dsdt_minust_center_68[i], 1000), quadrupole_form(np.linspace(0.01, phi_d_2H_dsdt_minust_center_68[i], 1000), curve_fit_params[0], curve_fit_params[1]), '-',  color = color_list[i-3])
    plt.plot(np.linspace(phi_d_2H_dsdt_minust_center_68[i], 0.6, 1000),  quadrupole_form(np.linspace(phi_d_2H_dsdt_minust_center_68[i], 0.6, 1000), curve_fit_params[0],  curve_fit_params[1]), ':', color = color_list[i-3])

    plt.axes(axs[1])
    this_radius = 0.1973*np.sqrt(24/curve_fit_params[1])
    this_radius_err = 0.1973*np.sqrt(24)*0.5*np.sqrt(curve_fit_cov[1][1])/curve_fit_params[1]**(3/2)
    plt.errorbar(phi_d_2H_dsdt_minust_center_68[i], this_radius, yerr=this_radius_err, fmt='o', color = color_list[i-3])

plt.axes(axs[0])
plt.yscale('log')
plt.xlabel(r'$-t\ [GeV^2/c]$')
plt.ylabel(r'$d\sigma/dt\ [nb/(GeV^2/c)]$')
plt.xlim(0, 0.6)
plt.ylim(1e0, 1e3)
plt.legend()

plt.axes(axs[1])
plt.xlabel(r'$-t\ [GeV^2/c]$')
plt.ylabel(r'$r\ [fm]$')
plt.xlim(0, 0.6)
plt.ylim(0, 3.0)

file_pdf.savefig()
plt.close()

# GAUSSIAN FORM
fig = plt.figure(figsize=(16, 6), dpi=300)
gs = fig.add_gridspec(1, 2)
axs = gs.subplots()
twin_axes = axs[1].twinx()
plt.suptitle("Gaussian Form, 5.8-7.8 GeV")
plt.axes(axs[0])
plt.errorbar(phi_d_2H_dsdt_minust_center_68[1:], phi_d_2H_dsdt_results_68[1:], yerr=phi_d_2H_dsdt_results_statserr_68[1:], fmt='b.', label='Data, 5.8-7.8 GeV')

for i in range(6, 16):
    curve_fit_params, curve_fit_cov = curve_fit(gaussian_form, \
                                            phi_d_2H_dsdt_minust_center_68[1:i], \
                                            phi_d_2H_dsdt_results_68[1:i], \
                                            sigma=phi_d_2H_dsdt_results_statserr_68[1:i], \
                                            absolute_sigma=True, p0=[300, 0.05], maxfev=50000)
    curve_fit_residuals             = phi_d_2H_dsdt_results_68[1:i] - gaussian_form(phi_d_2H_dsdt_minust_center_68[1:i], curve_fit_params[0], curve_fit_params[1])
    reduced_chi2                    = np.sum((curve_fit_residuals/phi_d_2H_dsdt_results_statserr_68[1:i])**2)/(len(phi_d_2H_dsdt_results_68[1:i])-2)
    plt.axes(axs[0])
    plt.plot(np.linspace(0.01, phi_d_2H_dsdt_minust_center_68[i], 1000), gaussian_form(np.linspace(0.01, phi_d_2H_dsdt_minust_center_68[i], 1000), curve_fit_params[0], curve_fit_params[1]), '-',  color = color_list[i-3])
    plt.plot(np.linspace(phi_d_2H_dsdt_minust_center_68[i], 0.6, 1000),  gaussian_form(np.linspace(phi_d_2H_dsdt_minust_center_68[i], 0.6, 1000), curve_fit_params[0],  curve_fit_params[1]), ':', color = color_list[i-3])

    plt.axes(axs[1])
    this_radius     = 0.1973*np.sqrt(3/2/curve_fit_params[1])
    this_radius_err = 0.1973*np.sqrt(3/2)*0.5*np.sqrt(curve_fit_cov[1][1])/curve_fit_params[1]**(3/2)
    if (i == 6):
        plt.errorbar(phi_d_2H_dsdt_minust_center_68[i], this_radius, yerr=this_radius_err, fmt='o', color = 'k', label='Extracted radius')
    else:
        plt.errorbar(phi_d_2H_dsdt_minust_center_68[i], this_radius, yerr=this_radius_err, fmt='o', color = 'k')

    plt.axes(twin_axes)
    if (i == 6):
        plt.errorbar(phi_d_2H_dsdt_minust_center_68[i], reduced_chi2, yerr=0, fmt='o', color = 'r', label='χ²/NDF')
    else:
        plt.plot(phi_d_2H_dsdt_minust_center_68[i], reduced_chi2, 'o', color = 'r')

plt.axes(axs[0])
plt.yscale('log')
plt.xlabel(r'$-t\ [GeV^2/c]$')
plt.ylabel(r'$d\sigma/dt\ [nb/(GeV^2/c)]$')
plt.xlim(0, 0.6)
plt.ylim(1e0, 1e3)
plt.legend()

plt.axes(axs[1])
plt.xlabel(r'$-t\ [GeV^2/c]$')
plt.ylabel(r'$r\ [fm]$')
plt.xlim(0, 0.6)
plt.ylim(0, 3.0)
plt.legend(loc='upper left')

plt.axes(twin_axes)
plt.ylabel(r'$\chi^2/ndf$')
plt.yticks(np.arange(0, 3.5, 0.5))
plt.legend(loc='upper right')

file_pdf.savefig()
plt.close()

# curve_fit_params, curve_fit_cov = curve_fit(uniform_form, \
#                                             phi_d_2H_dsdt_minust_center[fit_indices], \
#                                             phi_d_2H_dsdt_results[fit_indices], \
#                                             sigma=phi_d_2H_dsdt_results_statserr[fit_indices], \
#                                             absolute_sigma=True, p0=[300, 0.05], maxfev=5000)
# curve_fit_residuals             = phi_d_2H_dsdt_results[fit_indices] - uniform_form(phi_d_2H_dsdt_minust_center[fit_indices], curve_fit_params[0], curve_fit_params[1])
# reduced_chi2                    = np.sum((curve_fit_residuals/phi_d_2H_dsdt_results_statserr[fit_indices])**2)/(len(phi_d_2H_dsdt_results[fit_indices])-2)
# plt.plot(np.linspace(0.01, 0.6, 100), uniform_form(np.linspace(0.01, 0.6, 100), curve_fit_params[0], curve_fit_params[1]), '--', color = 'purple', label=f'Uniform fit, χ²/ndf = {reduced_chi2:.2f}, r={0.1973*np.sqrt(3*curve_fit_params[1]*curve_fit_params[1]/5):.2f} fm')

# # Format the plot
# plt.xlabel(r'$-t\ [GeV^2/c]$')
# plt.ylabel(r'$d\sigma/dt\ [nb/(GeV^2/c)]$')
# plt.xlim(0, 0.6)
# plt.ylim(1e0, 1e3)
# plt.yscale('log')
# plt.legend()
# file_pdf.savefig()
# plt.close()

file_pdf.close()