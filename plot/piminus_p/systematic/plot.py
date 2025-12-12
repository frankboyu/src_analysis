import numpy as np
import matplotlib.pyplot as plt

energy_low = 6
energy_hgh = 8

# bin_edges = np.loadtxt('flux_corr_90208.txt')
bin_edges = np.loadtxt('lumi_summed_2H.txt')
tagh_flux_error = np.loadtxt('tagh_flux_syst.txt')
tagm_flux_error = np.loadtxt('tagm_flux_syst.txt')

for i in range(len(bin_edges[:,0])):
    if bin_edges[i,0] <= 125 or bin_edges[i,0] >= 228:
        bin_edges[i,4] = tagh_flux_error[int(bin_edges[i,1])-1,2]
    elif bin_edges[i,0] >= 126 and bin_edges[i,0] <= 227:
        bin_edges[i,4] = tagm_flux_error[int(bin_edges[i,1])-1,2]



# plt.plot(bin_edges[:,3], bin_edges[:,6]/bin_edges[:,5]*100, marker='o', linestyle='none')
plt.plot(bin_edges[0:126,3], bin_edges[0:126,4], color='orange')
plt.plot(bin_edges[126:228,3], bin_edges[126:228,4], color='purple')
plt.plot(bin_edges[228:,3], bin_edges[228:,4], color='cyan')
# plt.plot(bin_edges[:,3], bin_edges[:,4], color='black')
plt.fill_between([6, 8], 0, 10, color='red', alpha=0.3, label='Energy bin 1')
plt.fill_between([8, 9], 0, 10, color='blue', alpha=0.3, label='Energy bin 2')
plt.fill_between([9, 11], 0, 10, color='green', alpha=0.3, label='Energy bin 3')
# plt.legend()
plt.xlabel("Photon energy [GeV]")
plt.ylabel("PS normalization uncertainty [%]")
plt.ylim(0, 10)
plt.savefig("flux_syst.png")

# file_output = open("output/flux_syst.txt", "w")

# for i in range(125):
    # file_output.write('{:>3.0f}    {:>3.0f}    {:>13.10f}    {:>13.10f}    {:>13.10f}    {:>17.16e}    {:>17.16e}\n'.format(bin_edges[i, 0], bin_edges[i, 1], bin_edges[i, 2], bin_edges[i, 3], bin_edges[i, 4], bin_edges[i, 5], bin_edges[i, 6]))