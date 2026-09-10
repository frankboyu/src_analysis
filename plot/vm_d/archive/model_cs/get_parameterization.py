import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages

file_pdf = PdfPages("/work/halld2/home/boyu/src_analysis/plot/vm_d/output/plots_vm_d_parameterization.pdf")

mass_photon = 0.0  # GeV/c^2
mass_proton = 0.9382720813  # GeV/c^2
mass_phi = 1.019461  # GeV/c^2
s_min = (mass_proton + mass_phi)**2
energy_min = (s_min - mass_proton**2) / (2*mass_proton)

def kallen(s, m1, m2):
    return s**2 + m1**4 + m2**4 - 2*s*m1**2 - 2*s*m2**2 - 2*m1**2*m2**2

def fit_function(x, c):
    return c*x

def slope_function(s, alpha, beta):
    return 2*alpha*np.log(s) + beta

# dsdt
data = np.loadtxt('input/parameterization_dsdt.txt', delimiter=',')
error = np.loadtxt('input/parameterization_dsdt_sigma.txt', delimiter=',')
energy = data[:, 0]
s = mass_proton**2 + 2*mass_proton*energy
q_photon = np.sqrt(kallen(s, mass_proton**2, mass_photon**2))/(2*np.sqrt(s))
q_phi = np.sqrt(kallen(s, mass_proton**2, mass_phi**2))/(2*np.sqrt(s))
x = (q_phi/q_photon)**2
dsdt = data[:, 1]
dsdt_error = error[:, 1] - data[:, 1]

fit_mask = energy > 0
x_fit = x[fit_mask]
dsdt_fit = dsdt[fit_mask]
dsdt_error_fit = dsdt_error[fit_mask]

popt, pcov = curve_fit(fit_function, x_fit, dsdt_fit, sigma=dsdt_error_fit, absolute_sigma=True)
coefficients = popt
fit_y = fit_function(x_fit, *coefficients)
chisquared_per_dof = np.sum(((dsdt_fit - fit_y) / dsdt_error_fit)**2) / (len(dsdt_fit) - len(coefficients))

plt.figure(figsize=(6, 6), dpi=300)
energy_plot = np.linspace(1.5, 10.0, 1000)
s_plot = mass_proton**2 + 2*mass_proton*energy_plot
q_photon_plot = np.sqrt(kallen(s_plot, mass_proton**2, mass_photon**2))/(2*np.sqrt(s_plot))
q_phi_plot = np.sqrt(kallen(s_plot, mass_proton**2, mass_phi**2))/(2*np.sqrt(s_plot))
x_plot = (q_phi_plot/q_photon_plot)**2
dsdt_plot = fit_function(x_plot, *coefficients)
# dsdt_plot = fit_function(x_plot, 2.93)
plt.errorbar(energy[0], dsdt[0], yerr=dsdt_error[0], color='k', label='Bonn (1974)', fmt='s', capsize=5)
plt.errorbar(energy[1:4], dsdt[1:4], yerr=dsdt_error[1:4], color='k', label='SLAC (1973)', fmt='o', capsize=5)
plt.errorbar(energy[4:7], dsdt[4:7], yerr=dsdt_error[4:7], color='k', label='Daresbury (1982)', fmt='d', capsize=5)
plt.errorbar(energy[7:14], dsdt[7:14], yerr=dsdt_error[7:14], color='k', label='DESY (1978)', fmt='^', capsize=5)
plt.errorbar(energy[14], dsdt[14], yerr=dsdt_error[14], color='k', label='Cornell (1972)', fmt='v', capsize=5)
plt.plot(energy_plot, dsdt_plot, color='r', linestyle='--', label='Fitted Curve')
plt.title('Parameterization Fit')
plt.xlabel('Photon Energy (GeV)')
plt.ylabel(r'$d\sigma/dt(\theta=0^{\circ})(\mu b/GeV^2)$')
plt.text(6, 1.5, r'$C=%.5f\pm%.5f$' % (coefficients[0], np.sqrt(pcov[0][0])), fontsize=10, color='r', ha='left', va='top')
plt.text(6, 1.0, r'$\chi^2/dof=%.2f$' % chisquared_per_dof, fontsize=10, color='r', ha='left', va='top')
plt.xlim(1.0, 10.0)
plt.ylim(0, 3.5)
# plt.yscale('log')
plt.legend()
plt.grid()
file_pdf.savefig()
plt.show()
plt.close()

energy_theory = np.array([2.0597, 3.0729, 6.8739, 8.3011, 9.6774])
s_theory = mass_proton**2 + 2*mass_proton*energy_theory
q_photon_theory = np.sqrt(kallen(s_theory, mass_proton**2, mass_photon**2))/(2*np.sqrt(s_theory))
q_phi_theory = np.sqrt(kallen(s_theory, mass_proton**2, mass_phi**2))/(2*np.sqrt(s_theory))
x_theory = (q_phi_theory/q_photon_theory)**2
dsdt_theory = fit_function(x_theory, *coefficients)
print("Theoretical ds/dt values at energies", energy_theory, "are:", dsdt_theory)

# slope factor
s           = np.loadtxt('input/parameterization_slope.txt')[:, 0]
slope       = np.loadtxt('input/parameterization_slope.txt')[:, 1]
slope_error = np.loadtxt('input/parameterization_slope.txt')[:, 2]
energy_slope = (s - mass_proton**2) / (2*mass_proton)
plt.figure(figsize=(6, 6), dpi=300)
plt.errorbar(energy_slope, slope, yerr=slope_error, color='b', label='Data Points', fmt='o')

popt, pcov = curve_fit(slope_function, s, slope, sigma=slope_error, absolute_sigma=True)
coefficients = popt
fit_y = slope_function(s, *coefficients)
chisquared_per_dof = np.sum(((slope - fit_y) / slope_error)**2) / (len(slope) - len(coefficients))

energy_plot = np.linspace(energy_min, 20.0, 1000)
s_plot = mass_proton**2 + 2*mass_proton*energy_plot
plt.plot(energy_plot, slope_function(s_plot, *coefficients), color='r', linestyle='--', label='Fitted Curve')
plt.text(2, 3, r'$\alpha ^{\prime} =%.5f\pm%.5f$' % (coefficients[0], np.sqrt(pcov[0][0])), fontsize=10, color='r', ha='left', va='top')
plt.text(2, 2.5, r'$b_0 =%.5f\pm%.5f$' % (coefficients[1], np.sqrt(pcov[1][1])), fontsize=10, color='r', ha='left', va='top')
plt.text(2, 2.0, r'$\chi^2/dof=%.2f$' % chisquared_per_dof, fontsize=10, color='r', ha='left', va='top')
plt.title('Energy Dependence of the Slope Factor')
plt.xlabel('Photon Energy (GeV)')
plt.ylabel(r'b (GeV$^{-2}$)')
plt.xlim(1.0, 15.0)
plt.ylim(0, 10)
plt.grid()
file_pdf.savefig()
plt.show()
plt.close()

b_theory = slope_function(s_theory, *coefficients)
print("Theoretical slope factor values at energies", energy_theory, "are:", b_theory)

file_pdf.close()