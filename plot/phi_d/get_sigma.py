import numpy as np
import matplotlib.pyplot as plt

# Load data
data_2H = np.loadtxt('output/yield_phi_d_recon_data_2H.txt')
data_4He = np.loadtxt('output/yield_phi_d_recon_data_4He.txt')
data_12C = np.loadtxt('output/yield_phi_d_recon_data_12C.txt')

plt.plot(data_2H[0:5,2], data_2H[0:5,4], 'o-', label='6.5 GeV')
plt.plot(data_2H[5:10,2], data_2H[5:10,4], 'o-', label='7.5 GeV')
plt.plot(data_2H[10:15,2], data_2H[10:15,4], 'o-', label='8.5 GeV')
plt.plot(data_2H[15:20,2], data_2H[15:20,4], 'o-', label='9.5 GeV')
plt.plot(data_2H[20:25,2], data_2H[20:25,4], 'o-', label='10.5 GeV')
plt.legend()
plt.xlabel(r'$\theta_{c.m.}$ (deg)')
plt.ylabel(r'Yield')
plt.title(r'$\phi D$ yield on 2H')

plt.savefig('output/yield_2H.png')
plt.close()

plt.plot(data_4He[0:5,2], data_4He[0:5,4], 'o-', label='6.5 GeV')
plt.plot(data_4He[5:10,2], data_4He[5:10,4], 'o-', label='7.5 GeV')
plt.plot(data_4He[10:15,2], data_4He[10:15,4], 'o-', label='8.5 GeV')
plt.plot(data_4He[15:20,2], data_4He[15:20,4], 'o-', label='9.5 GeV')
plt.plot(data_4He[20:25,2], data_4He[20:25,4], 'o-', label='10.5 GeV')
plt.legend()
plt.xlabel(r'$\theta_{c.m.}$ (deg)')
plt.ylabel(r'Yield')
plt.title(r'$\phi D$ yield on 4He')

plt.savefig('output/yield_4He.png')
plt.close()

plt.plot(data_12C[0:5,2], data_12C[0:5,4], 'o-', label='6.5 GeV')
plt.plot(data_12C[5:10,2], data_12C[5:10,4], 'o-', label='7.5 GeV')
plt.plot(data_12C[10:15,2], data_12C[10:15,4], 'o-', label='8.5 GeV')
plt.plot(data_12C[15:20,2], data_12C[15:20,4], 'o-', label='9.5 GeV')
plt.plot(data_12C[20:25,2], data_12C[20:25,4], 'o-', label='10.5 GeV')
plt.legend()
plt.xlabel(r'$\theta_{c.m.}$ (deg)')
plt.ylabel(r'Yield')
plt.title(r'$\phi D$ yield on 12C')

plt.savefig('output/yield_12C.png')
plt.close()