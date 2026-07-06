import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

clas_70 = np.loadtxt("70.dat")
clas_95 = np.loadtxt("95.dat")

one_sigma = np.loadtxt("1sig.dat", delimiter=',')
two_sigma = np.loadtxt("2sig.dat", delimiter=',')

plt.figure(figsize=(8, 6), dpi=300)
plt.plot(clas_70[:, 1], clas_70[:, 2], label='70% CL, CLAS', color='black', linestyle='solid', alpha=1.0)
plt.plot(clas_95[:, 1], clas_95[:, 2], label='95% CL, CLAS', color='black', linestyle='dashed', alpha=1.0)
plt.plot(one_sigma[:, 0]+5, one_sigma[:, 1], label='70% CL, SRC-CT', color='blue', linestyle='solid', alpha=1.0)
plt.plot(two_sigma[:, 0]+5, two_sigma[:, 1], label='95% CL, SRC-CT', color='blue', linestyle='dashed', alpha=1.0)
# plt.fill(one_sigma[:, 0]+5, one_sigma[:, 1], alpha=1.0, color='black')
# plt.fill(two_sigma[:, 0]+5, two_sigma[:, 1], alpha=0.5, color='black')
plt.xlabel(r'$\sigma_{\phi N}\ [mb]$', size=14)
plt.ylabel(r'$\rm b_{\phi N}\ [GeV^{-2}]$', size=14)
plt.xlim(0, 80)
plt.ylim(0, 20)
plt.xticks(np.arange(0, 85, 20), size=12)
plt.yticks(np.arange(0, 25, 5),  size=12)
plt.legend()
plt.savefig('confidence_contours.png', dpi=300)
plt.close()