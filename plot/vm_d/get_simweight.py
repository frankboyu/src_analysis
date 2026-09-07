import argparse
import glob
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, fsolve
from matplotlib.backends.backend_pdf import PdfPages

parser = argparse.ArgumentParser()
parser.add_argument('iteration', type=int, nargs='?', default=0,
                    help='Simulation-weight iteration to use as the starting iteration.')
args = parser.parse_args()
iteration = args.iteration

rad_to_deg = 180/np.pi

def dsdt_func(minust, a1, b1, a2, b2):
    return a1*np.exp(-b1*minust) + a2*np.exp(-b2*minust)

def bin_center(minust_low, minust_high, a1, b1, a2, b2):
    minust_center = np.zeros(len(minust_low), dtype=float)
    for i in range(len(minust_low)):
        this_minust_low = minust_low[i]
        this_minust_high = minust_high[i]
        f_average = (np.exp(-b1*this_minust_low) - np.exp(-b1*this_minust_high))*a1/b1/(this_minust_high - this_minust_low) + (np.exp(-b2*this_minust_low) - np.exp(-b2*this_minust_high))*a2/b2/(this_minust_high - this_minust_low)
        this_minust_center = fsolve(lambda minust: a1*np.exp(-b1*minust) + a2*np.exp(-b2*minust) - f_average, (this_minust_low + this_minust_high)/2)[0]
        minust_center[i] = this_minust_center
    return minust_center

def lumi(energy_min, energy_max, length):
    lumi_table = np.loadtxt('/work/halld2/home/boyu/src_analysis/flux/output/deuterium/lumi_summed_deuterium.txt')
    length_total = 29.5

    integrated_lumi = np.zeros(energy_min.shape, dtype=float)
    for i in range(len(energy_min)):
        for j in range(len(lumi_table)):
            if (lumi_table[j,3] > energy_min[i]) and (lumi_table[j,3] < energy_max[i]):
                integrated_lumi[i] += lumi_table[j][5]

    return integrated_lumi/length_total * length

###################################################################### DATA YIELD #####################################################################################

# Read the bin edges
phi_d_2H_dsdt_energy_center         = np.loadtxt('output/yield_dsdt/yield_phi_d_exc_recon_data_ver12_dsdt_nominal.txt')[:,0]
phi_d_2H_dsdt_energy_width          = np.loadtxt('output/yield_dsdt/yield_phi_d_exc_recon_data_ver12_dsdt_nominal.txt')[:,1]
phi_d_2H_dsdt_energy_low            = np.loadtxt('output/yield_dsdt/yield_phi_d_exc_recon_data_ver12_dsdt_nominal.txt')[:,2]
phi_d_2H_dsdt_energy_high           = np.loadtxt('output/yield_dsdt/yield_phi_d_exc_recon_data_ver12_dsdt_nominal.txt')[:,3]
phi_d_2H_dsdt_minust_center         = np.loadtxt('output/yield_dsdt/yield_phi_d_exc_recon_data_ver12_dsdt_nominal.txt')[:,4]
phi_d_2H_dsdt_minust_width          = np.loadtxt('output/yield_dsdt/yield_phi_d_exc_recon_data_ver12_dsdt_nominal.txt')[:,5]
phi_d_2H_dsdt_minust_low            = np.loadtxt('output/yield_dsdt/yield_phi_d_exc_recon_data_ver12_dsdt_nominal.txt')[:,6]
phi_d_2H_dsdt_minust_high           = np.loadtxt('output/yield_dsdt/yield_phi_d_exc_recon_data_ver12_dsdt_nominal.txt')[:,7]
phi_d_2H_dsdt_yield_data            = np.loadtxt('output/yield_dsdt/yield_phi_d_exc_recon_data_ver12_dsdt_nominal.txt')[:,8]
phi_d_2H_dsdt_yield_data_statserr   = np.loadtxt('output/yield_dsdt/yield_phi_d_exc_recon_data_ver12_dsdt_nominal.txt')[:,9]

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

###################################################################### PERFORM ITERATION #####################################################################################

# Calculate the bin centers
if iteration == 0:
    phi_d_2H_dsdt_minust_center_iteration   = (phi_d_2H_dsdt_minust_low + phi_d_2H_dsdt_minust_high)/2
else:
    last_paras = np.loadtxt(f'output/table_simweight_iter{iteration-1}.txt')
    phi_d_2H_dsdt_minust_center_iteration   = bin_center(phi_d_2H_dsdt_minust_low, phi_d_2H_dsdt_minust_high, last_paras[0], last_paras[1], last_paras[2], last_paras[3])

# Simulation yield numbers
phi_d_2H_dsdt_yield_sim_iteration               = np.loadtxt(f'output/yield_dsdt/yield_phi_d_exc_recon_sim_ver12_07_dsdt_simweight_iter{iteration}.txt')[:,8]
phi_d_2H_dsdt_yield_sim_statser_iteration       = np.loadtxt(f'output/yield_dsdt/yield_phi_d_exc_recon_sim_ver12_07_dsdt_simweight_iter{iteration}.txt')[:,9]
phi_d_2H_dsdt_yield_tagged_iteration            = np.loadtxt(f'output/yield_dsdt/yield_phi_d_exc_thrown_tagged_ver12_07_dsdt_simweight_iter{iteration}.txt')[:,8]
phi_d_2H_dsdt_yield_tagged_statserr_iteration   = np.loadtxt(f'output/yield_dsdt/yield_phi_d_exc_thrown_tagged_ver12_07_dsdt_simweight_iter{iteration}.txt')[:,9]

# Calculate the efficiency and differential cross section
phi_d_2H_dsdt_efficiency_iteration              = phi_d_2H_dsdt_yield_sim_iteration/phi_d_2H_dsdt_yield_tagged_iteration
phi_d_2H_dsdt_efficiency_statserr_iteration     = phi_d_2H_dsdt_efficiency_iteration*np.sqrt((phi_d_2H_dsdt_yield_sim_statser_iteration/phi_d_2H_dsdt_yield_sim_iteration)**2 + (phi_d_2H_dsdt_yield_tagged_statserr_iteration/phi_d_2H_dsdt_yield_tagged_iteration)**2)
phi_d_2H_dsdt_results_iteration                 = phi_d_2H_dsdt_yield_data/phi_d_2H_dsdt_efficiency_iteration/lumi(phi_d_2H_dsdt_energy_low, phi_d_2H_dsdt_energy_high, 28)/(phi_d_2H_dsdt_minust_high-phi_d_2H_dsdt_minust_low)/0.489/1000
phi_d_2H_dsdt_results_statserr_iteration        = phi_d_2H_dsdt_results_iteration*np.sqrt((phi_d_2H_dsdt_yield_data_statserr/phi_d_2H_dsdt_yield_data)**2 + (phi_d_2H_dsdt_efficiency_statserr_iteration/phi_d_2H_dsdt_efficiency_iteration)**2)

# Plot the cross section
fig = plt.figure(figsize=(8, 6), dpi=300)
color_code = ['b', 'k', 'r']
plt.errorbar(phi_d_2H_dsdt_minust_center_iteration[index[0]:index[1]],       phi_d_2H_dsdt_results_iteration[index[0]:index[1]],          xerr=phi_d_2H_dsdt_minust_width[index[0]:index[1]],        yerr=phi_d_2H_dsdt_results_statserr_iteration[index[0]:index[1]],            fmt='b.', label='6-8 GeV')
plt.errorbar(phi_d_2H_dsdt_minust_center_iteration[index[1]:index[2]],       phi_d_2H_dsdt_results_iteration[index[1]:index[2]],          xerr=phi_d_2H_dsdt_minust_width[index[1]:index[2]],        yerr=phi_d_2H_dsdt_results_statserr_iteration[index[1]:index[2]],            fmt='k.', label='8-9 GeV')
plt.errorbar(phi_d_2H_dsdt_minust_center_iteration[index[2]:index[3]],       phi_d_2H_dsdt_results_iteration[index[2]:index[3]],          xerr=phi_d_2H_dsdt_minust_width[index[2]:index[3]],        yerr=phi_d_2H_dsdt_results_statserr_iteration[index[2]:index[3]],            fmt='r.', label='9-11 GeV')

# Fit the cross section with the function
fit_indices = np.where(phi_d_2H_dsdt_minust_center_iteration > 0.24)[0]
curve_fit_params, curve_fit_cov = curve_fit(dsdt_func, \
                                            phi_d_2H_dsdt_minust_center_iteration[fit_indices], \
                                            phi_d_2H_dsdt_results_iteration[fit_indices], \
                                            sigma=phi_d_2H_dsdt_results_statserr_iteration[fit_indices], \
                                            absolute_sigma=True, p0=[3000, 15, 15, 3])
curve_fit_residuals             = phi_d_2H_dsdt_results_iteration[fit_indices] - dsdt_func(phi_d_2H_dsdt_minust_center_iteration[fit_indices], curve_fit_params[0], curve_fit_params[1], curve_fit_params[2], curve_fit_params[3])
reduced_chi2                    = np.sum((curve_fit_residuals/phi_d_2H_dsdt_results_statserr_iteration[fit_indices])**2)/(len(phi_d_2H_dsdt_results_iteration[fit_indices])-4)
plt.plot(np.linspace(0, 2, 100), dsdt_func(np.linspace(0, 2, 100), curve_fit_params[0], curve_fit_params[1], curve_fit_params[2], curve_fit_params[3]), '--', color = 'g', label='Combined fit')
plt.text(0.01, 1, r'$a_1=%.5f\pm%.5f$' % (curve_fit_params[0], np.sqrt(curve_fit_cov[0][0])), fontsize=10, color='b', ha='left', va='top')
plt.text(0.01, 0.7, r'$b_1=%.5f\pm%.5f$' % (curve_fit_params[1], np.sqrt(curve_fit_cov[1][1])), fontsize=10, color='b', ha='left', va='top')
plt.text(0.01, 0.5, r'$a_2=%.5f\pm%.5f$' % (curve_fit_params[2], np.sqrt(curve_fit_cov[2][2])), fontsize=10, color='b', ha='left', va='top')
plt.text(0.01, 0.35, r'$b_2=%.5f\pm%.5f$' % (curve_fit_params[3], np.sqrt(curve_fit_cov[3][3])), fontsize=10, color='b', ha='left', va='top')
plt.text(0.01, 0.25, r'$\chi^2/ndf=%.2f$' % (reduced_chi2), fontsize=10, color='b', ha='left', va='top')
np.savetxt(f'output/table_simweight_iter{iteration}.txt', curve_fit_params)

# Format the plot
plt.title("Simulation weight iteration %d" % iteration)
plt.xlabel(r'$-t\ [GeV^2/c]$')
plt.ylabel(r'$d\sigma/dt\ [nb/(GeV^2/c)]$')
plt.xlim(0, 2)
plt.ylim(1e-1, 1e3)
plt.yscale('log')
plt.legend()
file_pdf = PdfPages("/work/halld2/home/boyu/src_analysis/plot/vm_d/output/plots_simweight_iteration%d.pdf" % iteration)
file_pdf.savefig()
file_pdf.close()
plt.savefig("/work/halld2/home/boyu/src_analysis/plot/vm_d/output/plots_simweight_iteration%d.png" % iteration,
            dpi=300, bbox_inches='tight')
plt.close()

# Check for convergence
tolerance = 1e-4
if iteration > 0:
    if np.max((curve_fit_params-last_paras)/last_paras) < tolerance:
        print("Converged!")

        # Combine the iteration plots only after convergence is reached.
        iteration_plot_files = sorted(
            glob.glob('/work/halld2/home/boyu/src_analysis/plot/vm_d/output/plots_simweight_iteration*.png'),
            key=lambda filename: int(filename.rsplit('iteration', 1)[1].rsplit('.png', 1)[0]))
        with PdfPages('/work/halld2/home/boyu/src_analysis/plot/vm_d/output/plots_simweight_iterations.pdf') as iteration_plots_pdf:
            for iteration_plot_file in iteration_plot_files:
                plot_image = plt.imread(iteration_plot_file)
                iteration_fig = plt.figure(figsize=(8, 6), dpi=300)
                iteration_ax = iteration_fig.add_axes([0, 0, 1, 1])
                iteration_ax.imshow(plot_image)
                iteration_ax.axis('off')
                iteration_plots_pdf.savefig(iteration_fig, bbox_inches='tight', pad_inches=0)
                plt.close(iteration_fig)
        for iteration_plot_file in iteration_plot_files:
            os.remove(iteration_plot_file)
        iteration_pdf_files = glob.glob(
            '/work/halld2/home/boyu/src_analysis/plot/vm_d/output/'
            'plots_simweight_iteration*.pdf'
        )
        for iteration_pdf_file in iteration_pdf_files:
            if iteration_pdf_file != '/work/halld2/home/boyu/src_analysis/plot/vm_d/output/plots_simweight_iterations.pdf':
                os.remove(iteration_pdf_file)

        paras_list = []
        for i in range(iteration+1):
            paras_list.append(np.loadtxt(f'output/table_simweight_iter{i}.txt'))
            os.remove(f'output/table_simweight_iter{i}.txt')
        paras_list = np.array(paras_list)
        np.savetxt(f'output/table_simweight_iterations.txt', paras_list)
        fig = plt.figure(figsize=(8, 6), dpi=300)
        plt.plot(np.arange(1, iteration+1), (paras_list[1:,0]-paras_list[0:iteration,0])/paras_list[0:iteration,0], 'o-', label='a1')
        plt.plot(np.arange(1, iteration+1), (paras_list[1:,1]-paras_list[0:iteration,1])/paras_list[0:iteration,1], 'o-', label='b1')
        plt.plot(np.arange(1, iteration+1), (paras_list[1:,2]-paras_list[0:iteration,2])/paras_list[0:iteration,2], 'o-', label='a2')
        plt.plot(np.arange(1, iteration+1), (paras_list[1:,3]-paras_list[0:iteration,3])/paras_list[0:iteration,3], 'o-', label='b2')
        plt.axhline(y=tolerance, color='r', linestyle='--')
        plt.xlabel('Iteration')
        plt.ylabel('Parameter value')
        plt.xticks(np.arange(len(paras_list)-1)+1)
        plt.yscale('symlog', linthresh=tolerance)
        plt.ylim(0, 1)
        plt.legend()
        file_pdf = PdfPages("/work/halld2/home/boyu/src_analysis/plot/vm_d/output/plots_simweight_convergence.pdf")
        file_pdf.savefig()
        file_pdf.close()
        plt.close()
