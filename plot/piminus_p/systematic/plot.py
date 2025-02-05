import numpy as np
import matplotlib.pyplot as plt

bin_edges = np.loadtxt('input/flux_corr_90208.txt')
tagh_flux_error = np.loadtxt('input/tagh_flux_syst.txt')
tagm_flux_error = np.loadtxt('input/tagm_flux_syst.txt')

plt.plot(tagh_flux_error[:,1], tagh_flux_error[:,2])
plt.xlabel("Photon energy [GeV]")
plt.ylabel("PS normalization uncertainty [%]")
plt.savefig("output/flux_syst.png")

# file_output = open("output/flux_syst.txt", "w")

# for i in range(125):
    # file_output.write('{:>3.0f}    {:>3.0f}    {:>13.10f}    {:>13.10f}    {:>13.10f}    {:>17.16e}    {:>17.16e}\n'.format(bin_edges[i, 0], bin_edges[i, 1], bin_edges[i, 2], bin_edges[i, 3], bin_edges[i, 4], bin_edges[i, 5], bin_edges[i, 6]))