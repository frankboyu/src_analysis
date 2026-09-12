import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, fsolve
from matplotlib.backends.backend_pdf import PdfPages

def bin_center(minust_low, minust_high, a1, b1, a2, b2):
    f_average = (np.exp(-b1*minust_low) - np.exp(-b1*minust_high))*a1/b1/(minust_high - minust_low) + (np.exp(-b2*minust_low) - np.exp(-b2*minust_high))*a2/b2/(minust_high - minust_low)
    minust_center = fsolve(lambda minust: a1*np.exp(-b1*minust) + a2*np.exp(-b2*minust) - f_average, (minust_low + minust_high)/2)[0]
    return minust_center

integrated_lumi = np.loadtxt('/work/halld2/home/boyu/src_analysis/flux/output/deuterium/lumi_summed_deuterium.txt')
simweight_paras = np.loadtxt('/work/halld2/home/boyu/src_analysis/plot/vm_d/output/table_simweight_iterations.txt')
kinematic_bins  = np.loadtxt('/work/halld2/home/boyu/src_analysis/plot/vm_d/input/bins_phi_d_dsdt.txt')
bin_centers     = np.zeros((kinematic_bins.shape[0], 2))

for i in range(kinematic_bins.shape[0]):
    emin        = kinematic_bins[i, 0]
    emax        = kinematic_bins[i, 1]
    minustmin   = kinematic_bins[i, 2]
    minustmax   = kinematic_bins[i, 3]
    average_energy_numerator = np.sum(integrated_lumi[(integrated_lumi[:,3] >= emin) & (integrated_lumi[:,3] < emax), 3] * integrated_lumi[(integrated_lumi[:,3] >= emin) & (integrated_lumi[:,3] < emax), 5])
    average_energy_denominator = np.sum(integrated_lumi[(integrated_lumi[:,3] >= emin) & (integrated_lumi[:,3] < emax), 5])
    average_energy = average_energy_numerator / average_energy_denominator if average_energy_denominator != 0 else 0
    bin_centers[i, 0] = average_energy
    bin_centers[i, 1] = bin_center(minustmin, minustmax, simweight_paras[-1, 0], simweight_paras[-1, 1], simweight_paras[-1, 2], simweight_paras[-1, 3])

np.savetxt('/work/halld2/home/boyu/src_analysis/plot/vm_d/output/table_bin_centers.txt', bin_centers, fmt='%.6f')
