import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

mass_photon = 0.0  # GeV/c^2
mass_proton = 0.9382720813  # GeV/c^2
mass_phi = 1.019461  # GeV/c^2
s_min = (mass_proton + mass_phi)**2

def kallen(s, m1, m2):
    return s**2 + m1**4 + m2**4 - 2*s*m1**2 - 2*s*m2**2 - 2*m1**2*m2**2

# def fit_function(x, a, b, c, d, e, f):
    # return a + b*x + c*x**2 + d*x**3 + e*x**4 + f*x**5

# def fit_function(x, a, b, c, d, e):
    # return (a + b*x + c*x**2 + d*x**3 + e*x**4)*(x-1.573)

# def fit_function(x, a, b, c, d):
    # return (a + b*x + c*x**2 + d*x**3)*(x-1.573)

def fit_function(x, c):
    return c*x

# Load data from the text file
data = np.loadtxt('dsdt.txt', delimiter=',')
error = np.loadtxt('sigma.txt', delimiter=',')
energy = data[:, 0]
s = mass_proton**2 + 2*mass_proton*energy
q_photon = np.sqrt(kallen(s, mass_proton**2, mass_photon**2))/(2*np.sqrt(s))
q_phi = np.sqrt(kallen(s, mass_proton**2, mass_phi**2))/(2*np.sqrt(s))
x = (q_phi/q_photon)**2
dsdt = data[:, 1]
dsdt_error = error[:, 1] - data[:, 1]

fit_mask = energy > 2.7
x_fit = x[fit_mask]
dsdt_fit = dsdt[fit_mask]
dsdt_error_fit = dsdt_error[fit_mask]

# Fit the data using numpy's polyfit for a 5th degree polynomial
# popt, pcov = curve_fit(fit_function, x, y, p0=[-0.69937E+01, 0.82732E+01, -0.32595E+01, 0.67068E+00, -0.69670E-01, 0.28879E-02], sigma=error[:, 1]-data[:, 1], absolute_sigma=True)
# popt, pcov = curve_fit(fit_function, x, y, p0=[-0.69937E+01, 0.82732E+01, -0.32595E+01, 0.67068E+00, -0.69670E-01], sigma=error[:, 1]-data[:, 1])
# popt, pcov = curve_fit(fit_function, x, y, p0=[-0.69937E+01, 0.82732E+01, -0.32595E+01, 0.67068E+00], sigma=error[:, 1]-data[:, 1])
popt, pcov = curve_fit(fit_function, x_fit, dsdt_fit, sigma=dsdt_error_fit, absolute_sigma=True)
# coefficients = np.array([-0.69937E+01, 0.82732E+01, -0.32595E+01, 0.67068E+00, -0.69670E-01, 0.28879E-02])
coefficients = popt
fit_y = fit_function(x_fit, *coefficients)
chisquared_per_dof = np.sum(((dsdt_fit - fit_y) / dsdt_error_fit)**2) / (len(dsdt_fit) - len(coefficients))
print(f"Chi-squared per degree of freedom: {chisquared_per_dof:.2f}")
# print("Fitted coefficients:", coefficients)
print("Fitted coefficient:", coefficients[0])
print("Uncertainties:", np.sqrt(np.diag(pcov)))

# Create the plot
plt.figure(figsize=(8, 6), dpi=300)
energy_plot = np.linspace(1.5, 10.0, 1000)
s_plot = mass_proton**2 + 2*mass_proton*energy_plot
q_photon_plot = np.sqrt(kallen(s_plot, mass_proton**2, mass_photon**2))/(2*np.sqrt(s_plot))
q_phi_plot = np.sqrt(kallen(s_plot, mass_proton**2, mass_phi**2))/(2*np.sqrt(s_plot))
x_plot = (q_phi_plot/q_photon_plot)**2
dsdt_plot = fit_function(x_plot, *coefficients)
# dsdt_plot = fit_function(x_plot, 2.93)
plt.errorbar(energy, dsdt, yerr=dsdt_error, color='b', label='Data Points', fmt='o')
plt.plot(energy_plot, dsdt_plot, color='r', linestyle='--', label='Fitted Curve')
plt.title('Parameterization Fit')
plt.xlabel('Photon Energy (GeV)')
plt.ylabel(r'$d\sigma/dt(\theta=0^{\circ})(\mu b/GeV^2)$')
plt.text(2, 3, r'$C=%.5f\pm%.5f$' % (coefficients[0], np.sqrt(pcov[0][0])), fontsize=10, color='r', ha='left', va='top')
plt.text(2, 2.5, r'$\chi^2/dof=%.2f$' % chisquared_per_dof, fontsize=10, color='r', ha='left', va='top')
# plt.xlim(, 30.0)
# plt.ylim(0.1, 4.0)
# plt.yscale('log')
plt.legend()
plt.grid()
plt.savefig('parameterization_fit.png')
plt.show()