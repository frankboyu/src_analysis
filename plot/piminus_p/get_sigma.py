import numpy as np
import matplotlib.pyplot as plt

RadToDeg = 180.0/np.pi
mass_neutron = 0.93956542052
mass_proton  = 0.93827208816
mass_piminus = 0.13957039
mass_photon  = 0

data_2H_missprot = np.loadtxt('output/yield_piminus_p_recon_data_2H_missprot.txt')[:,4]
data_4He_misshe3 = np.loadtxt('output/yield_piminus_p_recon_data_4He_misshe3.txt')[:,4]
data_12C_missb11 = np.loadtxt('output/yield_piminus_p_recon_data_12C_missb11.txt')[:,4]

sim_2H_missprot_flat = np.loadtxt('output/yield_piminus_p_recon_sim_2H_missprot_flat.txt')[:,4]
sim_4He_misshe3_flat = np.loadtxt('output/yield_piminus_p_recon_sim_4He_misshe3_flat.txt')[:,4]
sim_12C_missb11_flat = np.loadtxt('output/yield_piminus_p_recon_sim_12C_missb11_flat.txt')[:,4]

gen_2H_flat = np.loadtxt('output/yield_piminus_p_thrown_gen_2H_flat.txt')[:,4]
gen_4He_flat = np.loadtxt('output/yield_piminus_p_thrown_gen_4He_flat.txt')[:,4]
gen_12C_flat = np.loadtxt('output/yield_piminus_p_thrown_gen_12C_flat.txt')[:,4]

tagged_2H_flat = np.loadtxt('output/yield_piminus_p_thrown_tagged_2H_flat.txt')[:,4]
tagged_4He_flat = np.loadtxt('output/yield_piminus_p_thrown_tagged_4He_flat.txt')[:,4]
tagged_12C_flat = np.loadtxt('output/yield_piminus_p_thrown_tagged_12C_flat.txt')[:,4]

transparency = 0.95

def dOmega(theta_low, theta_high):
    return 2*np.pi*(np.cos(theta_low/RadToDeg)-np.cos(theta_high/RadToDeg))

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

theta_low = np.loadtxt('output/bin_edges.txt')[:,0]
theta_high = np.loadtxt('output/bin_edges.txt')[:,1]
energy_low = np.loadtxt('output/bin_edges.txt')[:,2]
energy_high = np.loadtxt('output/bin_edges.txt')[:,3]
sqrts_low = np.sqrt(mass_neutron**2 + 2*mass_neutron*energy_low)
sqrts_high = np.sqrt(mass_neutron**2 + 2*mass_neutron*energy_high)

acceptance_2H_missprot = np.zeros(len(theta_low))
acceptance_4He_misshe3 = np.zeros(len(theta_low))
acceptance_12C_missb11 = np.zeros(len(theta_low))

dsdo_2H_missprot = np.zeros(len(theta_low))
dsdo_4He_misshe3 = np.zeros(len(theta_low))
dsdo_12C_missb11 = np.zeros(len(theta_low))

scaled_dsdt_2H_missprot = np.zeros(len(theta_low))
scaled_dsdt_4He_misshe3 = np.zeros(len(theta_low))
scaled_dsdt_12C_missb11 = np.zeros(len(theta_low))

for i in range(len(theta_low)):
    if (tagged_2H_flat[i] == 0) or (sim_2H_missprot_flat[i] == 0):
        continue
    acceptance_2H_missprot[i] = sim_2H_missprot_flat[i]/tagged_2H_flat[i]
    dsdo_2H_missprot[i] = data_2H_missprot[i]/acceptance_2H_missprot[i]/lumi(energy_low[i], energy_high[i], '2H')/dOmega(theta_low[i], theta_high[i])/transparency/1000000
    energy_com = mass_neutron**2 + 2*mass_neutron*(energy_low[i] + energy_high[i])/2
    pi_com = (energy_com*energy_com + mass_photon**4 + mass_neutron**4 - 2*mass_photon**2*mass_neutron**2 - 2*energy_com*mass_photon**2 - 2*energy_com*mass_neutron**2)/4/energy_com
    pi_com = np.sqrt(pi_com)
    pf_com = (energy_com*energy_com + mass_piminus**4 + mass_proton**4 - 2*mass_piminus**2*mass_proton**2 - 2*energy_com*mass_piminus**2 - 2*energy_com*mass_proton**2)/4/energy_com
    pf_com = np.sqrt(pf_com)
    scaled_dsdt_2H_missprot[i] = dsdo_2H_missprot[i] * energy_com**7 * np.pi/pi_com/pf_com * 1000 / 10000000
    
for i in range(len(theta_low)):
    if (tagged_4He_flat[i] == 0) or (sim_4He_misshe3_flat[i] == 0):
        continue
    acceptance_4He_misshe3[i] = sim_4He_misshe3_flat[i]/tagged_4He_flat[i]
    dsdo_4He_misshe3[i] = data_4He_misshe3[i]/acceptance_4He_misshe3[i]/lumi(energy_low[i], energy_high[i], '4He')/dOmega(theta_low[i], theta_high[i])/1000000
    energy_com = mass_neutron**2 + 2*mass_neutron*(energy_low[i] + energy_high[i])/2
    pi_com = (energy_com*energy_com + mass_photon**4 + mass_neutron**4 - 2*mass_photon**2*mass_neutron**2 - 2*energy_com*mass_photon**2 - 2*energy_com*mass_neutron**2)/4/energy_com
    pi_com = np.sqrt(pi_com)
    pf_com = (energy_com*energy_com + mass_piminus**4 + mass_proton**4 - 2*mass_piminus**2*mass_proton**2 - 2*energy_com*mass_piminus**2 - 2*energy_com*mass_proton**2)/4/energy_com
    pf_com = np.sqrt(pf_com)
    scaled_dsdt_4He_misshe3[i] = dsdo_4He_misshe3[i] * energy_com**7 * np.pi/pi_com/pf_com * 1000 / 10000000
    
for i in range(len(theta_low)):
    if (tagged_12C_flat[i] == 0) or (sim_12C_missb11_flat[i] == 0):
        continue
    acceptance_12C_missb11[i] = sim_12C_missb11_flat[i]/tagged_12C_flat[i]
    dsdo_12C_missb11[i] = data_12C_missb11[i]/acceptance_12C_missb11[i]/lumi(energy_low[i], energy_high[i], '12C')/dOmega(theta_low[i], theta_high[i])/1000000
    energy_com = mass_neutron**2 + 2*mass_neutron*(energy_low[i] + energy_high[i])/2
    pi_com = (energy_com*energy_com + mass_photon**4 + mass_neutron**4 - 2*mass_photon**2*mass_neutron**2 - 2*energy_com*mass_photon**2 - 2*energy_com*mass_neutron**2)/4/energy_com
    pi_com = np.sqrt(pi_com)
    pf_com = (energy_com*energy_com + mass_piminus**4 + mass_proton**4 - 2*mass_piminus**2*mass_proton**2 - 2*energy_com*mass_piminus**2 - 2*energy_com*mass_proton**2)/4/energy_com
    pf_com = np.sqrt(pf_com)
    scaled_dsdt_12C_missb11[i] = dsdo_12C_missb11[i] * energy_com**7 * np.pi/pi_com/pf_com * 1000 / 10000000
    
fig, ax = plt.subplots(1, 7, figsize=(40, 10), dpi=200)
for i in range(len(theta_low)):
    plot_id = int(theta_low[i]/5 - 2)
    if plot_id > 6:
        continue
    ax[plot_id].plot(sqrts_low[i], scaled_dsdt_2H_missprot[i], 'bo')
    
fig, ax = plt.subplots(1, 7, figsize=(40, 10), dpi=200)
for i in range(len(theta_low)):
    plot_id = int(theta_low[i]/5 - 2)
    if plot_id > 6:
        continue
    ax[plot_id].plot(sqrts_low[i], dsdo_4He_misshe3[i]/dsdo_2H_missprot[i]/2, 'bo')
for i in range(7):
    ax[i].set_ylim(0, 1)
    
edges_energy = np.loadtxt('input/edges_energy.txt')
edges_sqrts = np.sqrt(mass_neutron**2 + 2*mass_neutron*edges_energy)
edges_theta = np.loadtxt('input/edges_theta.txt')
bins_energy = len(edges_energy)-1
bins_sqrt_s = len(edges_sqrts)-1
bins_theta = len(edges_theta)-1

acceptance = np.zeros((bins_energy, bins_theta))
dsigma_dOmega = np.zeros((bins_energy, bins_theta))
scaled_dxs = np.zeros((bins_energy, bins_theta))

for i in range(bins_energy):
    for j in range(bins_theta):
        acceptance[i][j] = sim_2H_missprot_model[bins_theta*i+j][4]/tagged_2H_model[bins_theta*i+j][4]
        if (sim_2H_missprot_model[bins_theta*i+j][4] == 0):
            continue
        dsigma_dOmega[i][j] = data_2H_missprot[bins_theta*i+j][4]/lumi(edges_energy[i], edges_energy[i+1], '2H')/dOmega(edges_theta[j], edges_theta[j+1])/acceptance[i][j]/transparency/1000000
        energy_com = mass_neutron**2 + 2*mass_neutron*(edges_energy[i] + edges_energy[i+1])/2
        pi_com = (energy_com*energy_com + mass_photon**4 + mass_neutron**4 - 2*mass_photon**2*mass_neutron**2 - 2*energy_com*mass_photon**2 - 2*energy_com*mass_neutron**2)/4/energy_com
        pi_com = np.sqrt(pi_com)
        pf_com = (energy_com*energy_com + mass_piminus**4 + mass_proton**4 - 2*mass_piminus**2*mass_proton**2 - 2*energy_com*mass_piminus**2 - 2*energy_com*mass_proton**2)/4/energy_com
        pf_com = np.sqrt(pf_com)
        scaled_dxs[i][j] = dsigma_dOmega[i][j] * energy_com**7 * np.pi/pi_com/pf_com * 1000 / 10000000
        
plt.plot(edges_sqrts[:-1], scaled_dxs[:,4], label='data', marker='.')