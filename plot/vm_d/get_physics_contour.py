import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

rad_to_deg = 180/np.pi
file_pdf = PdfPages("/work/halld2/home/boyu/src_analysis/plot/vm_d/output/plots_vm_d_contour.pdf")

def legend_without_duplicate_labels(ax):
    handles, labels = ax.get_legend_handles_labels()
    unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
    ax.legend(*zip(*unique))

chisq_array = np.loadtxt("/work/halld2/home/boyu/src_analysis/plot/vm_d/output/table_vm_d_chisq_two_para_nominal.txt")
sphin_list = chisq_array[:, 0].reshape(100,60)[:,0]
bphin_list = chisq_array[:, 1].reshape(100,60)[0,:]
ndf=79
chi2_array = chisq_array[:, 2].reshape(100,60).T/ndf

fig = plt.figure(figsize=(8, 6), dpi=300)
plt.contourf(sphin_list, bphin_list, np.log(chi2_array), levels=1000, cmap='viridis')
cbar = plt.colorbar()
cbar.set_label(r'$\chi^2/NDF$')
plt.xlabel(r'$\sigma_{\phi N}\ [mb]$')
plt.ylabel(r'$\rm b_{\phi N}\ [GeV^{-2}]$')
plt.xlim(20, 40)
plt.ylim(8, 14)
file_pdf.savefig()
plt.close()

fig = plt.figure(figsize=(8, 6), dpi=300)
confidence_levels = [2.30, 6.18, 11.83]  # 68%, 95%, 99.7% for 2 parameters
best_fit_idx = np.unravel_index(np.argmin(chi2_array), chi2_array.shape)
print(sphin_list[best_fit_idx[1]], bphin_list[best_fit_idx[0]])
best_chi2 = chi2_array[best_fit_idx]
for level in confidence_levels:
    contour_level = best_chi2 + level/ndf
    plt.contour(sphin_list, bphin_list, chi2_array, levels=[contour_level], colors='black', linestyles='dashed')
    # plt.text(best_fit_idx[1]*5 + 5, best_fit_idx[0]*5 + 20 + level, f'{int(level*100)/100}', color='black')

variation_list = []
variation_list.append(["alpha_m0.05", "alpha_p0.05", "alpha_m0.10", "alpha_m0.15"])
variation_list.append(["t2_0.5", "t2_1.0", "t2_1.5"])
# variation_list.append(["gn_variation_10.0_3.0", "gn_variation_10.0_3.5", "gn_variation_10.0_4.0", "gn_variation_10.0_4.5", "gn_variation_10.0_5.0"])
# variation_list.append(["gn_variation_10.5_3.0", "gn_variation_10.5_3.5", "gn_variation_10.5_4.0", "gn_variation_10.5_4.5", "gn_variation_10.5_5.0"])
# variation_list.append(["gn_variation_11.0_3.0", "gn_variation_11.0_3.5", "gn_variation_11.0_4.0", "gn_variation_11.0_4.5", "gn_variation_11.0_5.0"])
# variation_list.append(["gn_variation_11.5_3.0", "gn_variation_11.5_3.5", "gn_variation_11.5_4.0", "gn_variation_11.5_4.5", "gn_variation_11.5_5.0"])
# variation_list.append(["gn_variation_12.0_3.0", "gn_variation_12.0_3.5", "gn_variation_12.0_4.0", "gn_variation_12.0_4.5", "gn_variation_12.0_5.0"])
variation_list.append(["gn_variation_10.0_3.0", "gn_variation_10.5_3.0", "gn_variation_11.0_3.0", "gn_variation_11.5_3.0", "gn_variation_12.0_3.0"])
variation_list.append(["gn_variation_10.0_3.5", "gn_variation_10.5_3.5", "gn_variation_11.0_3.5", "gn_variation_11.5_3.5", "gn_variation_12.0_3.5"])
variation_list.append(["gn_variation_10.0_4.0", "gn_variation_10.5_4.0", "gn_variation_11.0_4.0", "gn_variation_11.5_4.0", "gn_variation_12.0_4.0"])
variation_list.append(["gn_variation_10.0_4.5", "gn_variation_10.5_4.5", "gn_variation_11.0_4.5", "gn_variation_11.5_4.5", "gn_variation_12.0_4.5"])
variation_list.append(["gn_variation_10.0_5.0", "gn_variation_10.5_5.0", "gn_variation_11.0_5.0", "gn_variation_11.5_5.0", "gn_variation_12.0_5.0"])

color_list = ['blue', 'green', 'orange', 'purple', 'brown', 'pink', 'gray']

for i,variation in enumerate(variation_list):
    for var in variation:
        chisq_array_var = np.loadtxt(f"/work/halld2/home/boyu/src_analysis/plot/vm_d/output/table_vm_d_chisq_two_para_{var}.txt")
        if chisq_array_var.shape[0] == 2000:
            sphin_list_var = chisq_array_var[:, 0].reshape(50,40)[:,0]
            bphin_list_var = chisq_array_var[:, 1].reshape(50,40)[0,:]
            chi2_array_var = chisq_array_var[:, 2].reshape(50,40).T/ndf
        elif chisq_array_var.shape[0] == 6000:
            sphin_list_var = chisq_array_var[:, 0].reshape(100,60)[:,0]
            bphin_list_var = chisq_array_var[:, 1].reshape(100,60)[0,:]
            chi2_array_var = chisq_array_var[:, 2].reshape(100,60).T/ndf
        best_fit_idx_var = np.unravel_index(np.argmin(chi2_array_var), chi2_array_var.shape)
        plt.scatter(sphin_list_var[best_fit_idx_var[1]], bphin_list_var[best_fit_idx_var[0]], color=color_list[i], s=10, zorder=10)
        print(f"Variation: {var}, Best fit: sigma_phiN = {sphin_list_var[best_fit_idx_var[1]]}, b_phiN = {bphin_list_var[best_fit_idx_var[0]]}, Chi2/NDF = {chi2_array_var[best_fit_idx_var]}")

plt.xlim(20, 40)
plt.ylim(8, 14)
plt.xlabel(r'$\sigma_{\phi N}\ [mb]$')
plt.ylabel(r'$\rm b_{\phi N}\ [GeV^{-2}]$')
file_pdf.savefig()
plt.close()


fig = plt.figure(figsize=(8, 6), dpi=300)

clas_contour_coherent_16 = np.loadtxt("/work/halld2/home/boyu/src_analysis/plot/vm_d/input/data_clas_contour_coherent_16_26.dat")
clas_contour_coherent_26 = np.loadtxt("/work/halld2/home/boyu/src_analysis/plot/vm_d/input/data_clas_contour_coherent_26_36.dat")
clas_contour_incoherent_16 = np.loadtxt("/work/halld2/home/boyu/src_analysis/plot/vm_d/input/data_clas_contour_incoherent_16_26.dat", delimiter=',')
clas_contour_incoherent_26 = np.loadtxt("/work/halld2/home/boyu/src_analysis/plot/vm_d/input/data_clas_contour_incoherent_26_36.dat", delimiter=',')

confidence_levels = [2.30]  # 68%, 95%, 99.7% for 2 parameters
best_fit_idx = np.unravel_index(np.argmin(chi2_array), chi2_array.shape)
best_chi2 = chi2_array[best_fit_idx]

plt.xlim(0, 120)
plt.ylim(0, 30)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel(r'$\sigma_{\phi N}\ [\rm mb]$', fontsize=15)
plt.ylabel(r'$b_{\phi N}\ [\rm GeV^{-2}]$', fontsize=15)
plt.minorticks_on()
plt.tick_params(which='both', direction='in', top=True, right=True)
ax = plt.gca()
for spine in ax.spines.values():
    spine.set_linewidth(1.2)

for level in confidence_levels:
    contour_level = best_chi2 + level/ndf
    plt.contourf(sphin_list, bphin_list, chi2_array, levels=[best_chi2, contour_level], colors=['crimson'], zorder=5)
    plt.contour(sphin_list, bphin_list, chi2_array, levels=[contour_level], colors='darkred', linewidths=1.5, zorder=5)

alpha_value = 0.5

plt.axvspan(24, 52, color="orange", alpha=0.15, hatch="//", label='LEPS 1.5-2.4 GeV', zorder=1)

plt.fill([10, 12, 12, 10], [2, 2, 8, 8], facecolor='green', alpha=0.75, edgecolor='green', linestyle='-', lw=1.5, zorder=4)

plt.fill(clas_contour_coherent_16[:, 1],    clas_contour_coherent_16[:, 2],     facecolor='blue', alpha=0.50, edgecolor='blue', linestyle='-', lw=1.5, zorder=4)
plt.fill(clas_contour_incoherent_16[:, 0],  clas_contour_incoherent_16[:, 1],   facecolor='cyan', alpha=0.25, edgecolor='cyan', linestyle='-',  lw=1.5, zorder=3)


plt.fill(clas_contour_coherent_26[:, 1], clas_contour_coherent_26[:, 2], facecolor='purple', alpha=0.50, edgecolor='purple', linestyle='-', lw=1.5, zorder=4)
plt.fill(clas_contour_incoherent_26[:, 0], clas_contour_incoherent_26[:, 1], facecolor='magenta', alpha=0.25, edgecolor='magenta', linestyle='-', lw=1.5, zorder=3)


plt.text(35, 7.0, 'CLAS 1.6-2.6 GeV\n' + r"$\gamma d \rightarrow \phi d$", color='white', fontsize=10, zorder=4, rotation=40, ha='center')
plt.text(80, 17, 'CLAS 2.6-3.6 GeV\n' + r"$\gamma d \rightarrow \phi d$", color='white', fontsize=10, zorder=4, rotation=30, ha='center')
plt.text(50, 23.0, 'CLAS 1.6-2.6 GeV\n' + r"$\gamma d \rightarrow \phi p(n)$", color='black', fontsize=10, zorder=4, rotation=40, ha='center')
plt.text(100, 24.0, 'CLAS 2.6-3.6 GeV\n' + r"$\gamma d \rightarrow \phi p(n)$", color='black', fontsize=10, zorder=4, rotation=30, ha='center')
plt.text(37, 2.0, 'LEPS 1.5-2.4 GeV\n' + r"$\gamma A \rightarrow \phi X$", color='black', fontsize=10, zorder=4, ha='center')
plt.text(12, 0.3, "VMD, 1.6-10 GeV\n" + r"$\gamma p \rightarrow \phi p$", color='black', fontsize=10, zorder=4, ha='center')
plt.text(30, 10, 'This work', color='black', fontsize=10, zorder=4, rotation=45, ha='center')

file_pdf.savefig()
plt.close()

fig = plt.figure(figsize=(8, 6), dpi=300)

clas_contour_coherent_16 = np.loadtxt("/work/halld2/home/boyu/src_analysis/plot/vm_d/input/data_clas_contour_coherent_16_26.dat")
clas_contour_coherent_26 = np.loadtxt("/work/halld2/home/boyu/src_analysis/plot/vm_d/input/data_clas_contour_coherent_26_36.dat")
clas_contour_incoherent_16 = np.loadtxt("/work/halld2/home/boyu/src_analysis/plot/vm_d/input/data_clas_contour_incoherent_16_26.dat", delimiter=',')
clas_contour_incoherent_26 = np.loadtxt("/work/halld2/home/boyu/src_analysis/plot/vm_d/input/data_clas_contour_incoherent_26_36.dat", delimiter=',')

confidence_levels = [2.30]  # 68%, 95%, 99.7% for 2 parameters
best_fit_idx = np.unravel_index(np.argmin(chi2_array), chi2_array.shape)
best_chi2 = chi2_array[best_fit_idx]

plt.xlim(0, 120)
plt.ylim(0, 30)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel(r'$\sigma_{\phi N}\ [\rm mb]$', fontsize=15)
plt.ylabel(r'$b_{\phi N}\ [\rm GeV^{-2}]$', fontsize=15)
plt.minorticks_on()
plt.tick_params(which='both', direction='in', top=True, right=True)
ax = plt.gca()
for spine in ax.spines.values():
    spine.set_linewidth(1.2)

for level in confidence_levels:
    contour_level = best_chi2 + level/ndf
    plt.contourf(sphin_list, bphin_list, chi2_array, levels=[best_chi2, contour_level], colors=['crimson'], zorder=5)
    plt.contour(sphin_list, bphin_list, chi2_array, levels=[contour_level], colors='darkred', linewidths=1.5, zorder=5)
    plt.plot(1000,1000, color='red', linewidth=1.5, label='This work')

alpha_value = 0.5

plt.fill([10, 12, 12, 10], [2, 2, 8, 8], facecolor='green', alpha=0.75, edgecolor='green', linestyle='-', lw=1.5, zorder=4, label="VMD, 1.6-10 GeV, " + r"$\gamma p \rightarrow \phi p$")

plt.axvspan(24, 52, color="orange", alpha=0.15, hatch="//", label='LEPS 1.5-2.4 GeV, ' + r"$\gamma A \rightarrow \phi X$", zorder=1)

plt.fill(clas_contour_coherent_16[:, 1],    clas_contour_coherent_16[:, 2],     facecolor='blue', alpha=0.50, edgecolor='blue', linestyle='-', lw=1.5, zorder=4, label='CLAS 1.6-2.6 GeV, ' + r"$\gamma d \rightarrow \phi d$")
plt.fill(clas_contour_incoherent_16[:, 0],  clas_contour_incoherent_16[:, 1],   facecolor='cyan', alpha=0.25, edgecolor='cyan', linestyle='-',  lw=1.5, zorder=3, label='CLAS 1.6-2.6 GeV, ' + r"$\gamma d \rightarrow \phi p(n)$")


plt.fill(clas_contour_coherent_26[:, 1], clas_contour_coherent_26[:, 2], facecolor='purple', alpha=0.50, edgecolor='purple', linestyle='-', lw=1.5, zorder=4, label='CLAS 2.6-3.6 GeV, ' + r"$\gamma d \rightarrow \phi d$")
plt.fill(clas_contour_incoherent_26[:, 0], clas_contour_incoherent_26[:, 1], facecolor='magenta', alpha=0.25, edgecolor='magenta', linestyle='-', lw=1.5, zorder=3, label='CLAS 2.6-3.6 GeV, ' + r"$\gamma d \rightarrow \phi p(n)$")


# plt.text(35, 7.0, 'CLAS 1.6-2.6 GeV\n' + r"$\gamma d \rightarrow \phi d$", color='white', fontsize=10, zorder=4, rotation=40, ha='center')
# plt.text(80, 17, 'CLAS 2.6-3.6 GeV\n' + r"$\gamma d \rightarrow \phi d$", color='white', fontsize=10, zorder=4, rotation=30, ha='center')
# plt.text(50, 23.0, 'CLAS 1.6-2.6 GeV\n' + r"$\gamma d \rightarrow \phi p(n)$", color='black', fontsize=10, zorder=4, rotation=40, ha='center')
# plt.text(100, 24.0, 'CLAS 2.6-3.6 GeV\n' + r"$\gamma d \rightarrow \phi p(n)$", color='black', fontsize=10, zorder=4, rotation=30, ha='center')
# plt.text(37, 2.0, 'LEPS 1.5-2.4 GeV\n' + r"$\gamma A \rightarrow \phi X$", color='black', fontsize=10, zorder=4, ha='center')
# plt.text(12, 0.3, "VMD, 1.6-10 GeV\n" + r"$\gamma p \rightarrow \phi p$", color='black', fontsize=10, zorder=4, ha='center')
# plt.text(30, 10, 'This work', color='black', fontsize=10, zorder=4, rotation=45, ha='center')

plt.legend(loc='lower right', fontsize=10, frameon=True, edgecolor='black', framealpha=0.8)

file_pdf.savefig()
plt.close()







file_pdf.close()