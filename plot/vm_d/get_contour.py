import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages

rad_to_deg = 180/np.pi
file_pdf = PdfPages("/work/halld2/home/boyu/src_analysis/plot/vm_d/output/plots_vm_d_chisq.pdf")

def legend_without_duplicate_labels(ax):
    handles, labels = ax.get_legend_handles_labels()
    unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
    ax.legend(*zip(*unique))






fig = plt.figure(figsize=(8, 6), dpi=300)
plt.contourf(sphin_list, bphin_list, chi2_array, levels=50, cmap='viridis')
cbar = plt.colorbar()
cbar.set_label(r'$\chi^2/NDF$')
plt.xlabel(r'$\sigma_{\phi N}\ [mb]$')
plt.ylabel(r'$\rm b_{\phi N}\ [GeV^{-2}]$')
file_pdf.savefig()
plt.close()

fig = plt.figure(figsize=(8, 6), dpi=300)
confidence_levels = [2.30, 6.18, 11.83]  # 68%, 95%, 99.7% for 2 parameters
best_fit_idx = np.unravel_index(np.argmin(chi2_array), chi2_array.shape)
best_chi2 = chi2_array[best_fit_idx]
for level in confidence_levels:
    contour_level = best_chi2 + level/ndf
    plt.contour(sphin_list, bphin_list, chi2_array, levels=[contour_level], colors='black', linestyles='dashed')
    # plt.text(best_fit_idx[1]*5 + 5, best_fit_idx[0]*5 + 20 + level, f'{int(level*100)/100}', color='black')
plt.xlabel(r'$\sigma_{\phi N}\ [mb]$')
plt.ylabel(r'$\rm b_{\phi N}\ [GeV^{-2}]$')
file_pdf.savefig()
plt.close()






# # sphin_list = np.arange(0,120,1)
# # bphin_list = np.arange(0,30,1)

# # fig = plt.figure(figsize=(8, 6), dpi=300)
# # # plt.contourf(sphin_list, bphin_list, chi2_array, levels=50, cmap='viridis')
# # # cbar = plt.colorbar()
# # # cbar.set_label(r'$\chi^2/NDF$')
# # confidence_levels = [2.30, 6.18]  # 68%, 95%, 99.7% for 2 parameters
# # best_fit_idx = np.unravel_index(np.argmin(chi2_array), chi2_array.shape)
# # best_chi2 = chi2_array[best_fit_idx]
# # for level in confidence_levels:
# #     contour_level = best_chi2 + level/ndf
# #     plt.contour(sphin_list, bphin_list, chi2_array, levels=[contour_level], colors='black', linestyles='dashed')
# #     # plt.text(best_fit_idx[1]*5 + 5, best_fit_idx[0]*5 + 20 + level, f'{int(level*100)/100}', color='black')
# # plt.xlabel(r'$\sigma_{\phi N}\ [mb]$')
# # plt.ylabel(r'$\rm b_{\phi N}\ [GeV^{-2}]$')
# # plt.title(r'$\chi^2/NDF$ map for $\phi-D$ scattering parameters')
# # file_pdf.savefig()
# # plt.close()






file_pdf.close()