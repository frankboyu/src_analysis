import numpy as np
import matplotlib.pyplot as plt

# Load data
data_2H = np.loadtxt('output/yield_rho_d_recon_data_2H.txt')[:,4]
sim_2H = np.loadtxt('output/yield_rho_d_recon_sim_2H.txt')[:,4]
tagged_2H = np.loadtxt('output/yield_rho_d_thrown_tagged_2H.txt')[:,4]

def lumi(energy_low, energy_high, target):
    if target == '2H':
        lumi_table = np.loadtxt('/work/halld2/home/boyu/src_analysis/flux/output/2H/lumi_summed_2H.txt')
    elif target == '4He':
        lumi_table = np.loadtxt('/work/halld2/home/boyu/src_analysis/flux/output/4He/lumi_summed_4He.txt')
    elif target == '12C':
        lumi_table = np.loadtxt('/work/halld2/home/boyu/src_analysis/flux/output/12C/lumi_summed_12C.txt')

    integrated_lumi = 0
    for i in range(len(lumi_table)):
        if (lumi_table[i,3] > energy_low) and (lumi_table[i,3] < energy_high):
            integrated_lumi += lumi_table[i][5]

    return integrated_lumi

minust_low = np.loadtxt('output/bin_edges.txt')[:,0]
minust_high = np.loadtxt('output/bin_edges.txt')[:,1]
energy_low = np.loadtxt('output/bin_edges.txt')[:,2]
energy_high = np.loadtxt('output/bin_edges.txt')[:,3]

acceptance_2H = np.zeros(len(minust_low))
dsdt_2H = np.zeros(len(minust_low))
err_2H = 1/np.sqrt(data_2H)

for i in range(len(minust_low)):
    if (sim_2H[i] == 0) or (tagged_2H[i] == 0):
        continue
    acceptance_2H[i] = sim_2H[i]/tagged_2H[i]
    dsdt_2H[i] = data_2H[i]/acceptance_2H[i]/lumi(energy_low[i], energy_high[i], '2H')/(minust_high[i]-minust_low[i])/1000

err_2H = dsdt_2H*err_2H

anderson_6gev = np.array(
[[0.1547911547911548, 815.99652953849],
[0.20294840294840294, 407.74590714368907],
[0.29926289926289923, 29.344888759617085],
[0.40245700245700244, 11.27062634729409],
[0.5022113022113022, 7.153942654376493],
[0.6054054054054054, 5.768182080714631],
[0.7051597051597052, 4.328761281083057],
[0.8049140049140049, 2.7476476475318843],
[0.908108108108108, 1.873817422860385],
[1.0078624078624079, 1.218187912010116],
[1.1076167076167076, 0.7732353283868117],
[1.203931203931204, 0.47920294733971397]]
)
anderson_6gev[:,1] = anderson_6gev[:,1] / 100 * 1000

anderson_12gev = np.array(
[[0.1547911547911548, 66.16325808990872],
[0.20294840294840294, 20.48924587155778],
[0.29926289926289923, 2.8812299600625297],
[0.40245700245700244, 1.0549072709492509],
[0.5022113022113022, 0.6232221168523759],
[0.601965601965602, 0.46770017506294687],
[0.7051597051597052, 0.35098795090695634],
[0.8049140049140049, 0.21237842199968443],
[0.9046683046683046, 0.13480555607075764],
[1.0044226044226043, 0.09193354900897525],
[1.1076167076167076, 0.07066238825659604],
[1.203931203931204, 0.0397958387382604],
[1.4034398034398035, 0.02241233023614643]]
)
anderson_12gev[:,1] = anderson_12gev[:,1] / 10 * 1000

anderson_18gev = np.array(
[[0.1479115479115479, 6.656014903509869],
[0.20294840294840294, 1.918467937631489],
[0.29926289926289923, 0.20245624089682518],
[0.4058968058968059, 0.10116550820144546],
[0.5056511056511056, 0.05697466997066967],
[0.6054054054054054, 0.0397958387382604],
[0.7017199017199017, 0.027796717062106947],
[0.8014742014742015, 0.018508453401584444],
[0.9012285012285012, 0.012622237966080401],
[1.0044226044226043, 0.008011867645750875],
[1.1076167076167076, 0.0054638654988185395],
[1.2073710073710073, 0.0037262006200280487]]
)
anderson_18gev[:,1] = anderson_18gev[:,1] * 1000

for i in range(9):
    plt.errorbar((minust_low[11*i:11*(i+1)]+minust_high[11*i:11*(i+1)])/2, data_2H[11*i:11*(i+1)], xerr=(minust_high[11*i:11*(i+1)]-minust_low[11*i:11*(i+1)])/2, yerr=np.sqrt(data_2H[11*i:11*(i+1)]), label=str(energy_low[11*i])+'-'+str(energy_high[11*i])+' GeV')
plt.legend()
plt.title(r"Yield of $d(\gamma, \rho d')$")
plt.xlabel(r'$-t[GeV^2/c]$')
plt.ylabel('Yield')
plt.savefig('output/fig_yield_2H.png', dpi=300)
plt.close()

for i in range(9):
    plt.errorbar((minust_low[11*i:11*(i+1)]+minust_high[11*i:11*(i+1)])/2, acceptance_2H[11*i:11*(i+1)], xerr=(minust_high[11*i:11*(i+1)]-minust_low[11*i:11*(i+1)])/2, label=str(energy_low[11*i])+'-'+str(energy_high[11*i])+' GeV')
plt.legend()
plt.title(r"Acceptance of $d(\gamma, \rho d')$")
plt.xlabel(r'$-t[GeV^2/c]$')
plt.ylabel('Acceptance')
plt.savefig('output/fig_acceptance_2H.png', dpi=300)
plt.close()

for i in range(9):
    plt.errorbar((minust_low[11*i:11*(i+1)]+minust_high[11*i:11*(i+1)])/2, dsdt_2H[11*i:11*(i+1)], yerr=err_2H[11*i:11*(i+1)], fmt='*-', label=str(energy_low[11*i])+'-'+str(energy_high[11*i])+' GeV')
plt.legend()
plt.title(r"Differential cross section of $d(\gamma, \rho d')$")
plt.xlabel(r'$-t\:[GeV^2/c]$')
plt.ylabel(r'$d\sigma/dt\:[nb/(GeV^2/c)]$')
plt.yscale('log')
plt.savefig('output/fig_dsdt_2H_bare.png', dpi=300)
plt.close()


plt.errorbar((minust_low[0:11]+minust_high[0:11])/2, dsdt_2H[0:11], yerr=err_2H[0:11], fmt='*-', label='SRC-CT 6.0-6.5 GeV')
plt.errorbar((minust_low[88:99]+minust_high[88:99])/2, dsdt_2H[88:99], yerr=err_2H[88:99], fmt='*-', label='SRC-CT 10.0-10.5 GeV')
plt.plot(anderson_6gev[:,0], anderson_6gev[:,1], 'o-', label='Anderson 6 GeV')
plt.plot(anderson_12gev[:,0], anderson_12gev[:,1], 'o-', label='Anderson 12 GeV')
plt.plot(anderson_18gev[:,0], anderson_18gev[:,1], 'o-', label='Anderson 18 GeV')
plt.legend()
plt.title(r"Differential cross section of $d(\gamma, \rho d')$")
plt.xlabel(r'$-t\:[GeV^2/c]$')
plt.ylabel(r'$d\sigma/dt\:[nb/(GeV^2/c)]$')
plt.yscale('log')
plt.savefig('output/fig_dsdt_2H_compare.png', dpi=300)
plt.close()