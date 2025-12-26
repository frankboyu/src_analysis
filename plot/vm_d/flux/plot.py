import numpy as np
import matplotlib.pyplot as plt

energy_bins = [(6,8), (8,9), (9,11)]

region1 = 0
region2 = 125
region3 = 227

bin_edges = np.loadtxt('/work/halld2/home/boyu/src_analysis/flux/output/2H/lumi_summed_2H.txt')
tagh_flux_error = np.loadtxt('/work/halld2/home/boyu/src_analysis/flux/systematic/tagh_flux_syst.txt')
tagm_flux_error = np.loadtxt('/work/halld2/home/boyu/src_analysis/flux/systematic/tagm_flux_syst.txt')

for i in range(len(bin_edges[:,0])):
    if i < region2 or i >= region3:
        bin_edges[i,6] = tagh_flux_error[int(bin_edges[i,1])-1,2]
    elif i >= region2 and i < region3:
        bin_edges[i,6] = tagm_flux_error[int(bin_edges[i,1])-1,2]

fig, ax1 = plt.subplots(figsize=(10,6), dpi=300)

ax1.plot(bin_edges[region1:region2,3], bin_edges[region1:region2,6], color='orange', label='TAGH (upstream)')
ax1.plot(bin_edges[region2:region3,3], bin_edges[region2:region3,6], color='purple', label='TAGM')
ax1.plot(bin_edges[region3:,3], bin_edges[region3:,6], color='cyan', label='TAGH (downstream)')
ax1.fill_between([energy_bins[0][0], energy_bins[0][1]], 0, 10, color='red', alpha=0.3, label='Energy bin 1')
ax1.fill_between([energy_bins[1][0], energy_bins[1][1]], 0, 10, color='blue', alpha=0.3, label='Energy bin 2')
ax1.fill_between([energy_bins[2][0], energy_bins[2][1]], 0, 10, color='green', alpha=0.3, label='Energy bin 3')
ax1.legend(fontsize=14)
ax1.set_xlabel("Photon energy [GeV]", fontsize=16)
ax1.set_ylabel("PS normalization uncertainty [%]", fontsize=16)
ax1.set_ylim(0, 10)
ax1.set_xlim(6, 11)
ax1.tick_params(axis='x', labelsize=14)
ax1.tick_params(axis='y', labelsize=14)

ax2 = ax1.twinx()
ax2.errorbar(bin_edges[:,3], bin_edges[:,5], xerr=(bin_edges[:,4]-bin_edges[:,2])/2, fmt='.', color='black', label='Lumi per energy bin')
ax2.set_ylabel("Integrated luminosity [nb$^{-1}$]")


plt.savefig("flux_syst.png")

print(np.max(bin_edges[:,6]))

for i in range(3):
    energy_low, energy_high = energy_bins[i]
    flux_error = np.average(bin_edges[:,6][(bin_edges[:,3] >= energy_low) & (bin_edges[:,3] < energy_high)], weights=bin_edges[:,5][(bin_edges[:,3] >= energy_low) & (bin_edges[:,3] < energy_high)])
    print(f"Average flux uncertainty between {energy_low} and {energy_high} GeV: {flux_error:.2f} %")