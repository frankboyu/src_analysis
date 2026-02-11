import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages

rad_to_deg = 180/np.pi
file_pdf = PdfPages("/work/halld2/home/boyu/src_analysis/plot/vm_d/output/plots_vm_d.pdf")

def dsdt_func(minust, a1, b1, a2, b2):
    return a1*np.exp(-b1*minust) + a2*np.exp(-b2*minust)

def Wcostheta_func(costheta, c):
    return 0.75*((3*c-1)*costheta**2 + (1-c))

def Wphi_func(phi, c):
    return 1-2*c*np.cos(2*phi*np.pi/180)

def lumi(energy_min, energy_max, length):
    lumi_table = np.loadtxt('/work/halld2/home/boyu/src_analysis/flux/output/2H/lumi_summed_2H.txt')
    length_total = 29.5

    integrated_lumi = np.zeros(energy_min.shape, dtype=float)
    for i in range(len(energy_min)):
        for j in range(len(lumi_table)):
            if (lumi_table[j,3] > energy_min[i]) and (lumi_table[j,3] < energy_max[i]):
                integrated_lumi[i] += lumi_table[j][5]

    return integrated_lumi/length_total * length

def flux_systematic(energy_low, energy_high):
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

    flux_error = np.average(bin_edges[:,6][(bin_edges[:,3] >= energy_low) & (bin_edges[:,3] < energy_high)], weights=bin_edges[:,5][(bin_edges[:,3] >= energy_low) & (bin_edges[:,3] < energy_high)])
    return flux_error/100

def normalize_distribution(results, results_statserr, energy_bins, t_bins):
    for i in range(len(results)):
        if (i == 0):
            index = 0
        elif (i == len(results) - 1):
                N_total                     = np.sum(results[index:i+1])
                N_this                      = results[index:i+1]
                N_rest                      = np.sum(results[index:i+1]) - N_this
                sigma_this                  = results_statserr[index:i+1]**2
                sigma_rest                  = np.sqrt(np.sum(results_statserr[index:i+1]**2) - sigma_this**2)
                results[index:i+1]          = results[index:i+1]/N_total
                results_statserr[index:i+1] = np.sqrt(N_rest**2 * sigma_this**2 + N_this**2 * sigma_rest**2)/N_total**2
        else:
            if (energy_bins[i] != energy_bins[i-1]) or (t_bins[i] != t_bins[i-1]):
                N_total                     = np.sum(results[index:i])
                N_this                      = results[index:i]
                N_rest                      = np.sum(results[index:i]) - N_this
                sigma_this                  = results_statserr[index:i]**2
                sigma_rest                  = np.sqrt(np.sum(results_statserr[index:i]**2) - sigma_this**2)
                results[index:i]            = results[index:i]/N_total
                results_statserr[index:i]   = np.sqrt(N_rest**2 * sigma_this**2 + N_this**2 * sigma_rest**2)/N_total**2
                index = i
    return results, results_statserr

def legend_without_duplicate_labels(ax):
    handles, labels = ax.get_legend_handles_labels()
    unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
    ax.legend(*zip(*unique))

# List of labels and colors for the plots
dsdt_label_list = [r'$E_\gamma$=5.8-7.8 GeV', r'$E_\gamma$=7.8-8.8 GeV', r'$E_\gamma$=8.8-10.8 GeV']
dsdt_color_list = ['r', 'g', 'b']
W_label_list    = [ r'$E_\gamma$=5.8-7.8 GeV, $-t<0.6 \rm{ GeV}^2$', \
                    r'$E_\gamma$=5.8-7.8 GeV, $-t>0.6 \rm{ GeV}^2$', \
                    r'$E_\gamma$=7.8-8.8 GeV, $-t<0.6 \rm{ GeV}^2$', \
                    r'$E_\gamma$=7.8-8.8 GeV, $-t>0.6 \rm{ GeV}^2$', \
                    r'$E_\gamma$=8.8-10.8 GeV, $-t<0.6 \rm{ GeV}^2$', \
                    r'$E_\gamma$=8.8-10.8 GeV, $-t>0.6 \rm{ GeV}^2$']
W_color_list    = ['r', 'orange', 'green', 'cyan', 'blue', 'purple']
xlabel_list = [r'$-t[GeV^2/c]$', r'$\cos\vartheta$', r'$\varphi$[deg]', r'$\Phi$[deg]', r'$\psi$[deg]']
title_list  = [r'$d\sigma /dt$', r'$W(\cos\vartheta)$', r'$W(\varphi$)', r'$W(\Phi$)', r'$W(\psi$)']

color_options   = (['r', 'orange', 'black', 'cyan', 'blue', 'green', 'yellow', 'purple', 'brown', 'pink', 'gray', 'olive', 'navy'])
syst_list       = (['dEdx cut', 'Missing p- cut', 'Chi2/NDF cut', 'Momentum cut', 'Theta cut', 'Vertex Z cut', 'Vertex R cut', 'Beam accid subtraction', 'Combo accid subtraction', 'Fit max', 'Fit width', 'Fit background model', 'Fit signal model', 'Simulation weight'])
tag_list        = []
subplot_list    = []
color_list      = []

nominal_list    = (['nominal'])
dedx_list       = (['dEdx_1.50', 'dEdx_1.75', 'dEdx_2.50', 'dEdx_3.00'])
pminus_list     = (['misspminus_0.015', 'misspminus_0.0175', 'misspminus_0.025', 'misspminus_0.030'])
chisquared_list = (['chisquared_4.50', 'chisquared_4.75', 'chisquared_5.50', 'chisquared_6.00'])
momentum_list   = (['momentum_0.350', 'momentum_0.375', 'momentum_0.425', 'momentum_0.450'])
theta_list      = (['theta_1.90', 'theta_1.95', 'theta_2.05', 'theta_2.10'])
vertexZ_list    = (['vertexZ_13.50', 'vertexZ_13.75', 'vertexZ_14.25', 'vertexZ_14.50'])
vertexR_list    = (['vertexR_0.50', 'vertexR_0.75', 'vertexR_1.25', 'vertexR_1.50'])
beamaccid_list  = (['beamaccid_3', 'beamaccid_4out', 'beamaccid_5'])
comboaccid_list = (['comboaccid_all', 'comboaccid_none'])
fitmax_list     = (['fitmax_1.06', 'fitmax_1.07', 'fitmax_1.09', 'fitmax_1.10'])
fitwidth_list   = (['fitwidth_0.0040', 'fitwidth_0.0048', 'fitwidth_0.0060', 'fitwidth_0.0075'])
fitbkg_list     = (['fitbkg_fulllinear', 'fitbkg_quadratic', 'fitbkg_fullquadratic', 'fitbkg_phenomenological'])
fitsig_list     = (['fitsig_noBL', 'fitsig_nonrel'])

for i, this_list in enumerate([nominal_list, dedx_list, pminus_list, chisquared_list, momentum_list, theta_list, vertexZ_list, vertexR_list, beamaccid_list, comboaccid_list, fitmax_list, fitwidth_list, fitbkg_list, fitsig_list]):
    for j, this_label in enumerate(this_list):
            tag_list.append(this_label)
            subplot_list.append(i-1)
            color_list.append(color_options[j])

# Read the bin edges
phi_d_2H_dsdt_energy_center         = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_dsdt_nominal.txt')[:,0]
phi_d_2H_dsdt_energy_width          = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_dsdt_nominal.txt')[:,1]
phi_d_2H_dsdt_energy_low            = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_dsdt_nominal.txt')[:,2]
phi_d_2H_dsdt_energy_high           = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_dsdt_nominal.txt')[:,3]
phi_d_2H_dsdt_minust_center         = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_dsdt_nominal.txt')[:,4]
phi_d_2H_dsdt_minust_width          = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_dsdt_nominal.txt')[:,5]
phi_d_2H_dsdt_minust_low            = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_dsdt_nominal.txt')[:,6]
phi_d_2H_dsdt_minust_high           = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_dsdt_nominal.txt')[:,7]

phi_d_2H_Wcostheta_energy_center    = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wcostheta_nominal.txt')[:,0]
phi_d_2H_Wcostheta_energy_width     = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wcostheta_nominal.txt')[:,1]
phi_d_2H_Wcostheta_energy_low       = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wcostheta_nominal.txt')[:,2]
phi_d_2H_Wcostheta_energy_high      = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wcostheta_nominal.txt')[:,3]
phi_d_2H_Wcostheta_minust_center    = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wcostheta_nominal.txt')[:,4]
phi_d_2H_Wcostheta_minust_width     = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wcostheta_nominal.txt')[:,5]
phi_d_2H_Wcostheta_minust_low       = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wcostheta_nominal.txt')[:,6]
phi_d_2H_Wcostheta_minust_high      = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wcostheta_nominal.txt')[:,7]
phi_d_2H_Wcostheta_costheta_center  = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wcostheta_nominal.txt')[:,8]
phi_d_2H_Wcostheta_costheta_width   = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wcostheta_nominal.txt')[:,9]
phi_d_2H_Wcostheta_costheta_low     = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wcostheta_nominal.txt')[:,10]
phi_d_2H_Wcostheta_costheta_high    = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wcostheta_nominal.txt')[:,11]

phi_d_2H_Wdecayphi_energy_center    = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wdecayphi_nominal.txt')[:,0]
phi_d_2H_Wdecayphi_energy_width     = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wdecayphi_nominal.txt')[:,1]
phi_d_2H_Wdecayphi_energy_low       = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wdecayphi_nominal.txt')[:,2]
phi_d_2H_Wdecayphi_energy_high      = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wdecayphi_nominal.txt')[:,3]
phi_d_2H_Wdecayphi_minust_center    = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wdecayphi_nominal.txt')[:,4]
phi_d_2H_Wdecayphi_minust_width     = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wdecayphi_nominal.txt')[:,5]
phi_d_2H_Wdecayphi_minust_low       = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wdecayphi_nominal.txt')[:,6]
phi_d_2H_Wdecayphi_minust_high      = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wdecayphi_nominal.txt')[:,7]
phi_d_2H_Wdecayphi_decayphi_center  = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wdecayphi_nominal.txt')[:,8]
phi_d_2H_Wdecayphi_decayphi_width   = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wdecayphi_nominal.txt')[:,9]
phi_d_2H_Wdecayphi_decayphi_low     = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wdecayphi_nominal.txt')[:,10]
phi_d_2H_Wdecayphi_decayphi_high    = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wdecayphi_nominal.txt')[:,11]

phi_d_2H_Wpolphi_energy_center      = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpolphi_nominal.txt')[:,0]
phi_d_2H_Wpolphi_energy_width       = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpolphi_nominal.txt')[:,1]
phi_d_2H_Wpolphi_energy_low         = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpolphi_nominal.txt')[:,2]
phi_d_2H_Wpolphi_energy_high        = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpolphi_nominal.txt')[:,3]
phi_d_2H_Wpolphi_minust_center      = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpolphi_nominal.txt')[:,4]
phi_d_2H_Wpolphi_minust_width       = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpolphi_nominal.txt')[:,5]
phi_d_2H_Wpolphi_minust_low         = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpolphi_nominal.txt')[:,6]
phi_d_2H_Wpolphi_minust_high        = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpolphi_nominal.txt')[:,7]
phi_d_2H_Wpolphi_polphi_center      = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpolphi_nominal.txt')[:,8]
phi_d_2H_Wpolphi_polphi_width       = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpolphi_nominal.txt')[:,9]
phi_d_2H_Wpolphi_polphi_low         = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpolphi_nominal.txt')[:,10]
phi_d_2H_Wpolphi_polphi_high        = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpolphi_nominal.txt')[:,11]

phi_d_2H_Wpsi_energy_center         = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpsi_nominal.txt')[:,0]
phi_d_2H_Wpsi_energy_width          = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpsi_nominal.txt')[:,1]
phi_d_2H_Wpsi_energy_low            = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpsi_nominal.txt')[:,2]
phi_d_2H_Wpsi_energy_high           = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpsi_nominal.txt')[:,3]
phi_d_2H_Wpsi_minust_center         = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpsi_nominal.txt')[:,4]
phi_d_2H_Wpsi_minust_width          = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpsi_nominal.txt')[:,5]
phi_d_2H_Wpsi_minust_low            = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpsi_nominal.txt')[:,6]
phi_d_2H_Wpsi_minust_high           = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpsi_nominal.txt')[:,7]
phi_d_2H_Wpsi_psi_center            = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpsi_nominal.txt')[:,8]
phi_d_2H_Wpsi_psi_width             = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpsi_nominal.txt')[:,9]
phi_d_2H_Wpsi_psi_low               = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpsi_nominal.txt')[:,10]
phi_d_2H_Wpsi_psi_high              = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpsi_nominal.txt')[:,11]

# Read the yield numbers
phi_d_2H_dsdt_yield_data                    = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_dsdt_nominal.txt')[:,8]
phi_d_2H_dsdt_yield_data_statserr           = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_dsdt_nominal.txt')[:,9]
# phi_d_2H_dsdt_yield_sim                     = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_dsdt_nominal.txt')[:,8]
# phi_d_2H_dsdt_yield_sim_statserr            = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_dsdt_nominal.txt')[:,9]
# phi_d_2H_dsdt_yield_tagged                  = np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_dsdt_nominal.txt')[:,8]
# phi_d_2H_dsdt_yield_tagged_statserr         = np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_dsdt_nominal.txt')[:,9]
# phi_d_2H_dsdt_efficiency                    = phi_d_2H_dsdt_yield_sim/phi_d_2H_dsdt_yield_tagged
# phi_d_2H_dsdt_efficiency_statserr           = phi_d_2H_dsdt_efficiency*np.sqrt((phi_d_2H_dsdt_yield_sim_statserr/phi_d_2H_dsdt_yield_sim)**2 + (phi_d_2H_dsdt_yield_tagged_statserr/phi_d_2H_dsdt_yield_tagged)**2)
# phi_d_2H_dsdt_results                       = phi_d_2H_dsdt_yield_data/phi_d_2H_dsdt_efficiency/lumi(phi_d_2H_dsdt_energy_low, phi_d_2H_dsdt_energy_high, 28)/(phi_d_2H_dsdt_minust_high-phi_d_2H_dsdt_minust_low)/0.489/1000
# phi_d_2H_dsdt_results_statserr              = phi_d_2H_dsdt_results*np.sqrt((phi_d_2H_dsdt_yield_data_statserr/phi_d_2H_dsdt_yield_data)**2 + (phi_d_2H_dsdt_efficiency_statserr/phi_d_2H_dsdt_efficiency)**2)

phi_d_2H_Wcostheta_yield_data               = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wcostheta_nominal.txt')[:,12]
phi_d_2H_Wcostheta_yield_data_statserr      = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wcostheta_nominal.txt')[:,13]
# phi_d_2H_Wcostheta_yield_sim                = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_Wcostheta_nominal.txt')[:,12]
# phi_d_2H_Wcostheta_yield_sim_statserr       = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_Wcostheta_nominal.txt')[:,13]
# phi_d_2H_Wcostheta_yield_tagged             = np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_Wcostheta_nominal.txt')[:,12]
# phi_d_2H_Wcostheta_yield_tagged_statserr    = np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_Wcostheta_nominal.txt')[:,13]
# phi_d_2H_Wcostheta_efficiency               = phi_d_2H_Wcostheta_yield_sim/phi_d_2H_Wcostheta_yield_tagged
# phi_d_2H_Wcostheta_efficiency_statserr      = phi_d_2H_Wcostheta_efficiency*np.sqrt((phi_d_2H_Wcostheta_yield_sim_statserr/phi_d_2H_Wcostheta_yield_sim)**2 + (phi_d_2H_Wcostheta_yield_tagged_statserr/phi_d_2H_Wcostheta_yield_tagged)**2)
# phi_d_2H_Wcostheta_results                  = phi_d_2H_Wcostheta_yield_data/phi_d_2H_Wcostheta_efficiency  # raw results
# phi_d_2H_Wcostheta_results_statserr         = phi_d_2H_Wcostheta_results*np.sqrt((phi_d_2H_Wcostheta_yield_data_statserr/phi_d_2H_Wcostheta_yield_data)**2 + (phi_d_2H_Wcostheta_efficiency_statserr/phi_d_2H_Wcostheta_efficiency)**2)
# phi_d_2H_Wcostheta_results                  = normalize_distribution(phi_d_2H_Wcostheta_results, phi_d_2H_Wcostheta_results_statserr, phi_d_2H_Wcostheta_energy_low, phi_d_2H_Wcostheta_minust_low)[0] # normalize to have the sum equal to 1
# phi_d_2H_Wcostheta_results_statserr         = normalize_distribution(phi_d_2H_Wcostheta_results, phi_d_2H_Wcostheta_results_statserr, phi_d_2H_Wcostheta_energy_low, phi_d_2H_Wcostheta_minust_low)[1] # set the proper statistical uncertainties
# phi_d_2H_Wcostheta_results                  = phi_d_2H_Wcostheta_results/(phi_d_2H_Wcostheta_costheta_high - phi_d_2H_Wcostheta_costheta_low)  # normalize to have the integral equal to 1
# phi_d_2H_Wcostheta_results_statserr         = phi_d_2H_Wcostheta_results_statserr/(phi_d_2H_Wcostheta_costheta_high - phi_d_2H_Wcostheta_costheta_low)  # set the proper statistical uncertainties

phi_d_2H_Wdecayphi_yield_data               = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wdecayphi_nominal.txt')[:,12]
phi_d_2H_Wdecayphi_yield_data_statserr      = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wdecayphi_nominal.txt')[:,13]
# phi_d_2H_Wdecayphi_yield_sim                = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_Wdecayphi_nominal.txt')[:,12]
# phi_d_2H_Wdecayphi_yield_sim_statserr       = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_Wdecayphi_nominal.txt')[:,13]
# phi_d_2H_Wdecayphi_yield_tagged             = np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_Wdecayphi_nominal.txt')[:,12]
# phi_d_2H_Wdecayphi_yield_tagged_statserr    = np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_Wdecayphi_nominal.txt')[:,13]
# phi_d_2H_Wdecayphi_efficiency               = phi_d_2H_Wdecayphi_yield_sim/phi_d_2H_Wdecayphi_yield_tagged
# phi_d_2H_Wdecayphi_efficiency_statserr      = phi_d_2H_Wdecayphi_efficiency*np.sqrt((phi_d_2H_Wdecayphi_yield_sim_statserr/phi_d_2H_Wdecayphi_yield_sim)**2 + (phi_d_2H_Wdecayphi_yield_tagged_statserr/phi_d_2H_Wdecayphi_yield_tagged)**2)
# phi_d_2H_Wdecayphi_results                  = phi_d_2H_Wdecayphi_yield_data/phi_d_2H_Wdecayphi_efficiency  # raw results
# phi_d_2H_Wdecayphi_results_statserr         = phi_d_2H_Wdecayphi_results*np.sqrt((phi_d_2H_Wdecayphi_yield_data_statserr/phi_d_2H_Wdecayphi_yield_data)**2 + (phi_d_2H_Wdecayphi_efficiency_statserr/phi_d_2H_Wdecayphi_efficiency)**2)
# phi_d_2H_Wdecayphi_results                  = normalize_distribution(phi_d_2H_Wdecayphi_results, phi_d_2H_Wdecayphi_results_statserr, phi_d_2H_Wdecayphi_energy_low, phi_d_2H_Wdecayphi_minust_low)[0] # normalize to have the sum equal to 1
# phi_d_2H_Wdecayphi_results_statserr         = normalize_distribution(phi_d_2H_Wdecayphi_results, phi_d_2H_Wdecayphi_results_statserr, phi_d_2H_Wdecayphi_energy_low, phi_d_2H_Wdecayphi_minust_low)[1] # set the proper statistical uncertainties
# phi_d_2H_Wdecayphi_results                  = 2*np.pi*phi_d_2H_Wdecayphi_results/((phi_d_2H_Wdecayphi_decayphi_high - phi_d_2H_Wdecayphi_decayphi_low)/180*np.pi)  # normalize to have the integral equal to 2pi
# phi_d_2H_Wdecayphi_results_statserr         = 2*np.pi*phi_d_2H_Wdecayphi_results_statserr/((phi_d_2H_Wdecayphi_decayphi_high - phi_d_2H_Wdecayphi_decayphi_low)/180*np.pi)  # set the proper statistical uncertainties

phi_d_2H_Wpolphi_yield_data                 = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpolphi_nominal.txt')[:,12]
phi_d_2H_Wpolphi_yield_data_statserr        = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpolphi_nominal.txt')[:,13]
# phi_d_2H_Wpolphi_yield_sim                  = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_Wpolphi_nominal.txt')[:,12]
# phi_d_2H_Wpolphi_yield_sim_statserr         = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_Wpolphi_nominal.txt')[:,13]
# phi_d_2H_Wpolphi_yield_tagged               = np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_Wpolphi_nominal.txt')[:,12]
# phi_d_2H_Wpolphi_yield_tagged_statserr      = np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_Wpolphi_nominal.txt')[:,13]
# phi_d_2H_Wpolphi_efficiency                 = phi_d_2H_Wpolphi_yield_sim/phi_d_2H_Wpolphi_yield_tagged
# phi_d_2H_Wpolphi_efficiency_statserr        = phi_d_2H_Wpolphi_efficiency*np.sqrt((phi_d_2H_Wpolphi_yield_sim_statserr/phi_d_2H_Wpolphi_yield_sim)**2 + (phi_d_2H_Wpolphi_yield_tagged_statserr/phi_d_2H_Wpolphi_yield_tagged)**2)
# phi_d_2H_Wpolphi_results                    = phi_d_2H_Wpolphi_yield_data/phi_d_2H_Wpolphi_efficiency  # raw results
# phi_d_2H_Wpolphi_results_statserr           = phi_d_2H_Wpolphi_results*np.sqrt((phi_d_2H_Wpolphi_yield_data_statserr/phi_d_2H_Wpolphi_yield_data)**2 + (phi_d_2H_Wpolphi_efficiency_statserr/phi_d_2H_Wpolphi_efficiency)**2)
# phi_d_2H_Wpolphi_results                    = normalize_distribution(phi_d_2H_Wpolphi_results, phi_d_2H_Wpolphi_results_statserr, phi_d_2H_Wpolphi_energy_low, phi_d_2H_Wpolphi_minust_low)[0] # normalize to have the sum equal to 1
# phi_d_2H_Wpolphi_results_statserr           = normalize_distribution(phi_d_2H_Wpolphi_results, phi_d_2H_Wpolphi_results_statserr, phi_d_2H_Wpolphi_energy_low, phi_d_2H_Wpolphi_minust_low)[1] # set the proper statistical uncertainties
# phi_d_2H_Wpolphi_results                    = 2*np.pi*phi_d_2H_Wpolphi_results/((phi_d_2H_Wpolphi_polphi_high - phi_d_2H_Wpolphi_polphi_low)/180*np.pi)  # normalize to have the integral equal to 2pi
# phi_d_2H_Wpolphi_results_statserr           = 2*np.pi*phi_d_2H_Wpolphi_results_statserr/((phi_d_2H_Wpolphi_polphi_high - phi_d_2H_Wpolphi_polphi_low)/180*np.pi)  # set the proper statistical uncertainties

phi_d_2H_Wpsi_yield_data                    = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpsi_nominal.txt')[:,12]
phi_d_2H_Wpsi_yield_data_statserr           = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpsi_nominal.txt')[:,13]
# phi_d_2H_Wpsi_yield_sim                     = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_Wpsi_nominal.txt')[:,12]
# phi_d_2H_Wpsi_yield_sim_statserr            = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_Wpsi_nominal.txt')[:,13]
# phi_d_2H_Wpsi_yield_tagged                  = np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_Wpsi_nominal.txt')[:,12]
# phi_d_2H_Wpsi_yield_tagged_statserr         = np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_Wpsi_nominal.txt')[:,13]
# phi_d_2H_Wpsi_efficiency                    = phi_d_2H_Wpsi_yield_sim/phi_d_2H_Wpsi_yield_tagged
# phi_d_2H_Wpsi_efficiency_statserr           = phi_d_2H_Wpsi_efficiency*np.sqrt((phi_d_2H_Wpsi_yield_sim_statserr/phi_d_2H_Wpsi_yield_sim)**2 + (phi_d_2H_Wpsi_yield_tagged_statserr/phi_d_2H_Wpsi_yield_tagged)**2)
# phi_d_2H_Wpsi_results                       = phi_d_2H_Wpsi_yield_data/phi_d_2H_Wpsi_efficiency  # raw results
# phi_d_2H_Wpsi_results_statserr              = phi_d_2H_Wpsi_results*np.sqrt((phi_d_2H_Wpsi_yield_data_statserr/phi_d_2H_Wpsi_yield_data)**2 + (phi_d_2H_Wpsi_efficiency_statserr/phi_d_2H_Wpsi_efficiency)**2)
# phi_d_2H_Wpsi_results                       = normalize_distribution(phi_d_2H_Wpsi_results, phi_d_2H_Wpsi_results_statserr, phi_d_2H_Wpsi_energy_low, phi_d_2H_Wpsi_minust_low)[0] # normalize to have the sum equal to 1
# phi_d_2H_Wpsi_results_statserr              = normalize_distribution(phi_d_2H_Wpsi_results, phi_d_2H_Wpsi_results_statserr, phi_d_2H_Wpsi_energy_low, phi_d_2H_Wpsi_minust_low)[1] # set the proper statistical uncertainties
# phi_d_2H_Wpsi_results                       = 2*np.pi*phi_d_2H_Wpsi_results/((phi_d_2H_Wpsi_psi_high - phi_d_2H_Wpsi_psi_low)/180*np.pi)  # normalize to have the integral equal to 2pi
# phi_d_2H_Wpsi_results_statserr              = 2*np.pi*phi_d_2H_Wpsi_results_statserr/((phi_d_2H_Wpsi_psi_high - phi_d_2H_Wpsi_psi_low)/180*np.pi)  # set the proper statistical uncertainties

for tag in tag_list:
    if (tag == 'nominal'):
        continue
    if tag.find('simweight') != -1:
        phi_d_2H_dsdt_yield_data                = np.vstack((phi_d_2H_dsdt_yield_data,                      np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_dsdt_nominal.txt')[:,8]))
        phi_d_2H_dsdt_yield_data_statserr       = np.vstack((phi_d_2H_dsdt_yield_data_statserr,             np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_dsdt_nominal.txt')[:,9]))
        # phi_d_2H_Wcostheta_yield_data           = np.vstack((phi_d_2H_Wcostheta_yield_data,                 np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wcostheta_nominal.txt')[:,12]))
        # phi_d_2H_Wcostheta_yield_data_statserr  = np.vstack((phi_d_2H_Wcostheta_yield_data_statserr,        np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wcostheta_nominal.txt')[:,13]))
        # phi_d_2H_Wdecayphi_yield_data           = np.vstack((phi_d_2H_Wdecayphi_yield_data,                 np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wdecayphi_nominal.txt')[:,12]))
        # phi_d_2H_Wdecayphi_yield_data_statserr  = np.vstack((phi_d_2H_Wdecayphi_yield_data_statserr,        np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wdecayphi_nominal.txt')[:,13]))
        # phi_d_2H_Wpolphi_yield_data             = np.vstack((phi_d_2H_Wpolphi_yield_data,                   np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpolphi_nominal.txt')[:,12]))
        # phi_d_2H_Wpolphi_yield_data_statserr    = np.vstack((phi_d_2H_Wpolphi_yield_data_statserr,          np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpolphi_nominal.txt')[:,13]))
        # phi_d_2H_Wpsi_yield_data                = np.vstack((phi_d_2H_Wpsi_yield_data,                      np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpsi_nominal.txt')[:,12]))
        # phi_d_2H_Wpsi_yield_data_statserr       = np.vstack((phi_d_2H_Wpsi_yield_data_statserr,             np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpsi_nominal.txt')[:,13]))
    else:
        phi_d_2H_dsdt_yield_data                = np.vstack((phi_d_2H_dsdt_yield_data,                      np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_dsdt_'+tag+'.txt')[:,8]))
        phi_d_2H_dsdt_yield_data_statserr       = np.vstack((phi_d_2H_dsdt_yield_data_statserr,             np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_dsdt_'+tag+'.txt')[:,9]))
        # phi_d_2H_Wcostheta_yield_data           = np.vstack((phi_d_2H_Wcostheta_yield_data,                 np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wcostheta_'+tag+'.txt')[:,12]))
        # phi_d_2H_Wcostheta_yield_data_statserr  = np.vstack((phi_d_2H_Wcostheta_yield_data_statserr,        np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wcostheta_'+tag+'.txt')[:,13]))
        # phi_d_2H_Wdecayphi_yield_data           = np.vstack((phi_d_2H_Wdecayphi_yield_data,                 np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wdecayphi_'+tag+'.txt')[:,12]))
        # phi_d_2H_Wdecayphi_yield_data_statserr  = np.vstack((phi_d_2H_Wdecayphi_yield_data_statserr,        np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wdecayphi_'+tag+'.txt')[:,13]))
        # phi_d_2H_Wpolphi_yield_data             = np.vstack((phi_d_2H_Wpolphi_yield_data,                   np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpolphi_'+tag+'.txt')[:,12]))
        # phi_d_2H_Wpolphi_yield_data_statserr    = np.vstack((phi_d_2H_Wpolphi_yield_data_statserr,          np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpolphi_'+tag+'.txt')[:,13]))
        # phi_d_2H_Wpsi_yield_data                = np.vstack((phi_d_2H_Wpsi_yield_data,                      np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpsi_'+tag+'.txt')[:,12]))
        # phi_d_2H_Wpsi_yield_data_statserr       = np.vstack((phi_d_2H_Wpsi_yield_data_statserr,             np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpsi_'+tag+'.txt')[:,13]))
    # if tag.find('fit') != -1:
    #     phi_d_2H_dsdt_yield_sim                 = np.vstack((phi_d_2H_dsdt_yield_sim,                       np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_dsdt_nominal.txt')[:,8]))
    #     phi_d_2H_dsdt_yield_sim_statserr        = np.vstack((phi_d_2H_dsdt_yield_sim_statserr,              np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_dsdt_nominal.txt')[:,9]))
    #     phi_d_2H_Wcostheta_yield_sim            = np.vstack((phi_d_2H_Wcostheta_yield_sim,                  np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_Wcostheta_nominal.txt')[:,12]))
    #     phi_d_2H_Wcostheta_yield_sim_statserr   = np.vstack((phi_d_2H_Wcostheta_yield_sim_statserr,         np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_Wcostheta_nominal.txt')[:,13]))
    #     phi_d_2H_Wdecayphi_yield_sim            = np.vstack((phi_d_2H_Wdecayphi_yield_sim,                  np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_Wdecayphi_nominal.txt')[:,12]))
    #     phi_d_2H_Wdecayphi_yield_sim_statserr   = np.vstack((phi_d_2H_Wdecayphi_yield_sim_statserr,         np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_Wdecayphi_nominal.txt')[:,13]))
    #     phi_d_2H_Wpolphi_yield_sim              = np.vstack((phi_d_2H_Wpolphi_yield_sim,                    np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_Wpolphi_nominal.txt')[:,12]))
    #     phi_d_2H_Wpolphi_yield_sim_statserr     = np.vstack((phi_d_2H_Wpolphi_yield_sim_statserr,           np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_Wpolphi_nominal.txt')[:,13]))
    #     phi_d_2H_Wpsi_yield_sim                 = np.vstack((phi_d_2H_Wpsi_yield_sim,                       np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_Wpsi_nominal.txt')[:,12]))
    #     phi_d_2H_Wpsi_yield_sim_statserr        = np.vstack((phi_d_2H_Wpsi_yield_sim_statserr,              np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_Wpsi_nominal.txt')[:,13]))
    # else:
    #     phi_d_2H_dsdt_yield_sim                 = np.vstack((phi_d_2H_dsdt_yield_sim,                       np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_dsdt_'+tag+'.txt')[:,8]))
    #     phi_d_2H_dsdt_yield_sim_statserr        = np.vstack((phi_d_2H_dsdt_yield_sim_statserr,              np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_dsdt_'+tag+'.txt')[:,9]))
    #     phi_d_2H_Wcostheta_yield_sim            = np.vstack((phi_d_2H_Wcostheta_yield_sim,                  np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_Wcostheta_'+tag+'.txt')[:,12]))
    #     phi_d_2H_Wcostheta_yield_sim_statserr   = np.vstack((phi_d_2H_Wcostheta_yield_sim_statserr,         np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_Wcostheta_'+tag+'.txt')[:,13]))
    #     phi_d_2H_Wdecayphi_yield_sim            = np.vstack((phi_d_2H_Wdecayphi_yield_sim,                  np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_Wdecayphi_'+tag+'.txt')[:,12]))
    #     phi_d_2H_Wdecayphi_yield_sim_statserr   = np.vstack((phi_d_2H_Wdecayphi_yield_sim_statserr,         np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_Wdecayphi_'+tag+'.txt')[:,13]))
    #     phi_d_2H_Wpolphi_yield_sim              = np.vstack((phi_d_2H_Wpolphi_yield_sim,                    np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_Wpolphi_'+tag+'.txt')[:,12]))
    #     phi_d_2H_Wpolphi_yield_sim_statserr     = np.vstack((phi_d_2H_Wpolphi_yield_sim_statserr,           np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_Wpolphi_'+tag+'.txt')[:,13]))
    #     phi_d_2H_Wpsi_yield_sim                 = np.vstack((phi_d_2H_Wpsi_yield_sim,                       np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_Wpsi_'+tag+'.txt')[:,12]))
    #     phi_d_2H_Wpsi_yield_sim_statserr        = np.vstack((phi_d_2H_Wpsi_yield_sim_statserr,              np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_Wpsi_'+tag+'.txt')[:,13]))

    # phi_d_2H_dsdt_yield_tagged                  = np.vstack((phi_d_2H_dsdt_yield_tagged,                    np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_dsdt_nominal.txt')[:,8]))
    # phi_d_2H_dsdt_yield_tagged_statserr         = np.vstack((phi_d_2H_dsdt_yield_tagged_statserr,           np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_dsdt_nominal.txt')[:,9]))
    # phi_d_2H_Wcostheta_yield_tagged             = np.vstack((phi_d_2H_Wcostheta_yield_tagged,               np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_Wcostheta_nominal.txt')[:,12]))
    # phi_d_2H_Wcostheta_yield_tagged_statserr    = np.vstack((phi_d_2H_Wcostheta_yield_tagged_statserr,      np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_Wcostheta_nominal.txt')[:,13]))
    # phi_d_2H_Wdecayphi_yield_tagged             = np.vstack((phi_d_2H_Wdecayphi_yield_tagged,               np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_Wdecayphi_nominal.txt')[:,12]))
    # phi_d_2H_Wdecayphi_yield_tagged_statserr    = np.vstack((phi_d_2H_Wdecayphi_yield_tagged_statserr,      np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_Wdecayphi_nominal.txt')[:,13]))
    # phi_d_2H_Wpolphi_yield_tagged               = np.vstack((phi_d_2H_Wpolphi_yield_tagged,                 np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_Wpolphi_nominal.txt')[:,12]))
    # phi_d_2H_Wpolphi_yield_tagged_statserr      = np.vstack((phi_d_2H_Wpolphi_yield_tagged_statserr,        np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_Wpolphi_nominal.txt')[:,13]))
    # phi_d_2H_Wpsi_yield_tagged                  = np.vstack((phi_d_2H_Wpsi_yield_tagged,                    np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_Wpsi_nominal.txt')[:,12]))
    # phi_d_2H_Wpsi_yield_tagged_statserr         = np.vstack((phi_d_2H_Wpsi_yield_tagged_statserr,           np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_Wpsi_nominal.txt')[:,13]))

    # phi_d_2H_dsdt_efficiency                    = np.vstack((phi_d_2H_dsdt_efficiency,                      phi_d_2H_dsdt_yield_sim[-1,:]/phi_d_2H_dsdt_yield_tagged[-1,:]))
    # phi_d_2H_dsdt_efficiency_statserr           = np.vstack((phi_d_2H_dsdt_efficiency_statserr,             phi_d_2H_dsdt_efficiency[-1,:]*np.sqrt((phi_d_2H_dsdt_yield_sim_statserr[-1,:]/phi_d_2H_dsdt_yield_sim[-1,:])**2 + (phi_d_2H_dsdt_yield_tagged_statserr[-1,:]/phi_d_2H_dsdt_yield_tagged[-1,:])**2)))
    # phi_d_2H_Wcostheta_efficiency               = np.vstack((phi_d_2H_Wcostheta_efficiency,                 phi_d_2H_Wcostheta_yield_sim[-1,:]/phi_d_2H_Wcostheta_yield_tagged[-1,:]))
    # phi_d_2H_Wcostheta_efficiency_statserr      = np.vstack((phi_d_2H_Wcostheta_efficiency_statserr,        phi_d_2H_Wcostheta_efficiency[-1,:]*np.sqrt((phi_d_2H_Wcostheta_yield_sim_statserr[-1,:]/phi_d_2H_Wcostheta_yield_sim[-1,:])**2 + (phi_d_2H_Wcostheta_yield_tagged_statserr[-1,:]/phi_d_2H_Wcostheta_yield_tagged[-1,:])**2)))
    # phi_d_2H_Wdecayphi_efficiency               = np.vstack((phi_d_2H_Wdecayphi_efficiency,                 phi_d_2H_Wdecayphi_yield_sim[-1,:]/phi_d_2H_Wdecayphi_yield_tagged[-1,:]))
    # phi_d_2H_Wdecayphi_efficiency_statserr      = np.vstack((phi_d_2H_Wdecayphi_efficiency_statserr,        phi_d_2H_Wdecayphi_efficiency[-1,:]*np.sqrt((phi_d_2H_Wdecayphi_yield_sim_statserr[-1,:]/phi_d_2H_Wdecayphi_yield_sim[-1,:])**2 + (phi_d_2H_Wdecayphi_yield_tagged_statserr[-1,:]/phi_d_2H_Wdecayphi_yield_tagged[-1,:])**2)))
    # phi_d_2H_Wpolphi_efficiency                 = np.vstack((phi_d_2H_Wpolphi_efficiency,                   phi_d_2H_Wpolphi_yield_sim[-1,:]/phi_d_2H_Wpolphi_yield_tagged[-1,:]))
    # phi_d_2H_Wpolphi_efficiency_statserr        = np.vstack((phi_d_2H_Wpolphi_efficiency_statserr,          phi_d_2H_Wpolphi_efficiency[-1,:]*np.sqrt((phi_d_2H_Wpolphi_yield_sim_statserr[-1,:]/phi_d_2H_Wpolphi_yield_sim[-1,:])**2 + (phi_d_2H_Wpolphi_yield_tagged_statserr[-1,:]/phi_d_2H_Wpolphi_yield_tagged[-1,:])**2)))
    # phi_d_2H_Wpsi_efficiency                    = np.vstack((phi_d_2H_Wpsi_efficiency,                      phi_d_2H_Wpsi_yield_sim[-1,:]/phi_d_2H_Wpsi_yield_tagged[-1,:]))
    # phi_d_2H_Wpsi_efficiency_statserr           = np.vstack((phi_d_2H_Wpsi_efficiency_statserr,             phi_d_2H_Wpsi_efficiency[-1,:]*np.sqrt((phi_d_2H_Wpsi_yield_sim_statserr[-1,:]/phi_d_2H_Wpsi_yield_sim[-1,:])**2 + (phi_d_2H_Wpsi_yield_tagged_statserr[-1,:]/phi_d_2H_Wpsi_yield_tagged[-1,:])**2)))

    # phi_d_2H_dsdt_results                       = np.vstack((phi_d_2H_dsdt_results,                         phi_d_2H_dsdt_yield_data[-1,:]/phi_d_2H_dsdt_efficiency[-1,:]/lumi(phi_d_2H_dsdt_energy_low, phi_d_2H_dsdt_energy_high, 28)/(phi_d_2H_dsdt_minust_high-phi_d_2H_dsdt_minust_low)/0.489/1000))
    # phi_d_2H_dsdt_results_statserr              = np.vstack((phi_d_2H_dsdt_results_statserr,                phi_d_2H_dsdt_results[-1,:]*np.sqrt((phi_d_2H_dsdt_yield_data_statserr[-1,:]/phi_d_2H_dsdt_yield_data[-1,:])**2 + (phi_d_2H_dsdt_efficiency_statserr[-1,:]/phi_d_2H_dsdt_efficiency[-1,:])**2)))
    # phi_d_2H_Wcostheta_results                  = np.vstack((phi_d_2H_Wcostheta_results,                    phi_d_2H_Wcostheta_yield_data[-1,:]/phi_d_2H_Wcostheta_efficiency[-1,:]))
    # phi_d_2H_Wcostheta_results_statserr         = np.vstack((phi_d_2H_Wcostheta_results_statserr,           phi_d_2H_Wcostheta_results[-1,:]*np.sqrt((phi_d_2H_Wcostheta_yield_data_statserr[-1,:]/phi_d_2H_Wcostheta_yield_data[-1,:])**2 + (phi_d_2H_Wcostheta_efficiency_statserr[-1,:]/phi_d_2H_Wcostheta_efficiency[-1,:])**2)))
    # phi_d_2H_Wcostheta_results[-1,:]            = normalize_distribution(phi_d_2H_Wcostheta_results[-1,:],  phi_d_2H_Wcostheta_results_statserr[-1,:], phi_d_2H_Wcostheta_energy_low, phi_d_2H_Wcostheta_minust_low)[0]
    # phi_d_2H_Wcostheta_results_statserr[-1,:]   = normalize_distribution(phi_d_2H_Wcostheta_results[-1,:],  phi_d_2H_Wcostheta_results_statserr[-1,:], phi_d_2H_Wcostheta_energy_low, phi_d_2H_Wcostheta_minust_low)[1]
    # phi_d_2H_Wcostheta_results[-1,:]            = phi_d_2H_Wcostheta_results[-1,:]/                         (phi_d_2H_Wcostheta_costheta_high - phi_d_2H_Wcostheta_costheta_low)
    # phi_d_2H_Wcostheta_results_statserr[-1,:]   = phi_d_2H_Wcostheta_results_statserr[-1,:]/                (phi_d_2H_Wcostheta_costheta_high - phi_d_2H_Wcostheta_costheta_low)
    # phi_d_2H_Wdecayphi_results                  = np.vstack((phi_d_2H_Wdecayphi_results,                    phi_d_2H_Wdecayphi_yield_data[-1,:]/phi_d_2H_Wdecayphi_efficiency[-1,:]))
    # phi_d_2H_Wdecayphi_results_statserr         = np.vstack((phi_d_2H_Wdecayphi_results_statserr,           phi_d_2H_Wdecayphi_results[-1,:]*np.sqrt((phi_d_2H_Wdecayphi_yield_data_statserr[-1,:]/phi_d_2H_Wdecayphi_yield_data[-1,:])**2 + (phi_d_2H_Wdecayphi_efficiency_statserr[-1,:]/phi_d_2H_Wdecayphi_efficiency[-1,:])**2)))
    # phi_d_2H_Wdecayphi_results[-1,:]            = normalize_distribution(phi_d_2H_Wdecayphi_results[-1,:],  phi_d_2H_Wdecayphi_results_statserr[-1,:], phi_d_2H_Wdecayphi_energy_low, phi_d_2H_Wdecayphi_minust_low)[0]
    # phi_d_2H_Wdecayphi_results_statserr[-1,:]   = normalize_distribution(phi_d_2H_Wdecayphi_results[-1,:],  phi_d_2H_Wdecayphi_results_statserr[-1,:], phi_d_2H_Wdecayphi_energy_low, phi_d_2H_Wdecayphi_minust_low)[1]
    # phi_d_2H_Wdecayphi_results[-1,:]            = 2*np.pi*phi_d_2H_Wdecayphi_results[-1,:]/                 ((phi_d_2H_Wdecayphi_decayphi_high - phi_d_2H_Wdecayphi_decayphi_low)/180*np.pi)
    # phi_d_2H_Wdecayphi_results_statserr[-1,:]   = 2*np.pi*phi_d_2H_Wdecayphi_results_statserr[-1,:]/        ((phi_d_2H_Wdecayphi_decayphi_high - phi_d_2H_Wdecayphi_decayphi_low)/180*np.pi)
    # phi_d_2H_Wpolphi_results                    = np.vstack((phi_d_2H_Wpolphi_results,                      phi_d_2H_Wpolphi_yield_data[-1,:]/phi_d_2H_Wpolphi_efficiency[-1,:]))
    # phi_d_2H_Wpolphi_results_statserr           = np.vstack((phi_d_2H_Wpolphi_results_statserr,             phi_d_2H_Wpolphi_results[-1,:]*np.sqrt((phi_d_2H_Wpolphi_yield_data_statserr[-1,:]/phi_d_2H_Wpolphi_yield_data[-1,:])**2 + (phi_d_2H_Wpolphi_efficiency_statserr[-1,:]/phi_d_2H_Wpolphi_efficiency[-1,:])**2)))
    # phi_d_2H_Wpolphi_results[-1,:]              = normalize_distribution(phi_d_2H_Wpolphi_results[-1,:],    phi_d_2H_Wpolphi_results_statserr[-1,:], phi_d_2H_Wpolphi_energy_low, phi_d_2H_Wpolphi_minust_low)[0]
    # phi_d_2H_Wpolphi_results_statserr[-1,:]     = normalize_distribution(phi_d_2H_Wpolphi_results[-1,:],    phi_d_2H_Wpolphi_results_statserr[-1,:], phi_d_2H_Wpolphi_energy_low, phi_d_2H_Wpolphi_minust_low)[1]
    # phi_d_2H_Wpolphi_results[-1,:]              = 2*np.pi*phi_d_2H_Wpolphi_results[-1,:]/                   ((phi_d_2H_Wpolphi_polphi_high - phi_d_2H_Wpolphi_polphi_low)/180*np.pi)
    # phi_d_2H_Wpolphi_results_statserr[-1,:]     = 2*np.pi*phi_d_2H_Wpolphi_results_statserr[-1,:]/          ((phi_d_2H_Wpolphi_polphi_high - phi_d_2H_Wpolphi_polphi_low)/180*np.pi)
    # phi_d_2H_Wpsi_results                       = np.vstack((phi_d_2H_Wpsi_results,                         phi_d_2H_Wpsi_yield_data[-1,:]/phi_d_2H_Wpsi_efficiency[-1,:]))
    # phi_d_2H_Wpsi_results_statserr              = np.vstack((phi_d_2H_Wpsi_results_statserr,                phi_d_2H_Wpsi_results[-1,:]*np.sqrt((phi_d_2H_Wpsi_yield_data_statserr[-1,:]/phi_d_2H_Wpsi_yield_data[-1,:])**2 + (phi_d_2H_Wpsi_efficiency_statserr[-1,:]/phi_d_2H_Wpsi_efficiency[-1,:])**2)))
    # phi_d_2H_Wpsi_results[-1,:]                 = normalize_distribution(phi_d_2H_Wpsi_results[-1,:],       phi_d_2H_Wpsi_results_statserr[-1,:], phi_d_2H_Wpsi_energy_low, phi_d_2H_Wpsi_minust_low)[0]
    # phi_d_2H_Wpsi_results_statserr[-1,:]        = normalize_distribution(phi_d_2H_Wpsi_results[-1,:],       phi_d_2H_Wpsi_results_statserr[-1,:], phi_d_2H_Wpsi_energy_low, phi_d_2H_Wpsi_minust_low)[1]
    # phi_d_2H_Wpsi_results[-1,:]                 = 2*np.pi*phi_d_2H_Wpsi_results[-1,:]/                      ((phi_d_2H_Wpsi_psi_high - phi_d_2H_Wpsi_psi_low)/180*np.pi)

################################################################# INDICES OF THE KINEMATIC BINS ##########################################################################################################
dsdt_index = []
for i in range(len(phi_d_2H_dsdt_yield_data[0,:])):
    if (i == 0):
        dsdt_index.append(i)
    elif (i == len(phi_d_2H_dsdt_yield_data[0,:]) - 1):
        dsdt_index.append(i+1)
    else:
        if (phi_d_2H_dsdt_energy_low[i] != phi_d_2H_dsdt_energy_low[i-1]):
            dsdt_index.append(i)

# Wcostheta_index = []
# for i in range(len(phi_d_2H_Wcostheta_yield_data[0,:])):
#     if (i == 0):
#         Wcostheta_index.append(i)
#     elif (i == len(phi_d_2H_Wcostheta_yield_data[0,:]) - 1):
#         Wcostheta_index.append(i+1)
#     else:
#         if ((phi_d_2H_Wcostheta_energy_low[i] != phi_d_2H_Wcostheta_energy_low[i-1]) or (phi_d_2H_Wcostheta_minust_low[i] != phi_d_2H_Wcostheta_minust_low[i-1])):
#             Wcostheta_index.append(i)

# Wdecayphi_index = []
# for i in range(len(phi_d_2H_Wdecayphi_yield_data[0,:])):
#     if (i == 0):
#         Wdecayphi_index.append(i)
#     elif (i == len(phi_d_2H_Wdecayphi_yield_data[0,:]) - 1):
#         Wdecayphi_index.append(i+1)
#     else:
#         if ((phi_d_2H_Wdecayphi_energy_low[i] != phi_d_2H_Wdecayphi_energy_low[i-1]) or (phi_d_2H_Wdecayphi_minust_low[i] != phi_d_2H_Wdecayphi_minust_low[i-1])):
#             Wdecayphi_index.append(i)

# Wpolphi_index = []
# for i in range(len(phi_d_2H_Wpolphi_yield_data[0,:])):
#     if (i == 0):
#         Wpolphi_index.append(i)
#     elif (i == len(phi_d_2H_Wpolphi_yield_data[0,:]) - 1):
#         Wpolphi_index.append(i+1)
#     else:
#         if ((phi_d_2H_Wpolphi_energy_low[i] != phi_d_2H_Wpolphi_energy_low[i-1]) or (phi_d_2H_Wpolphi_minust_low[i] != phi_d_2H_Wpolphi_minust_low[i-1])):
#             Wpolphi_index.append(i)

# Wpsi_index = []
# for i in range(len(phi_d_2H_Wpsi_yield_data[0,:])):
#     if (i == 0):
#         Wpsi_index.append(i)
#     elif (i == len(phi_d_2H_Wpsi_yield_data[0,:]) - 1):
#         Wpsi_index.append(i+1)
#     else:
#         if ((phi_d_2H_Wpsi_energy_low[i] != phi_d_2H_Wpsi_energy_low[i-1]) or (phi_d_2H_Wpsi_minust_low[i] != phi_d_2H_Wpsi_minust_low[i-1])):
#             Wpsi_index.append(i)

################################################################# DATA YIELD ##############################################################################################################################
fig = plt.figure(figsize=(24, 16), dpi=300)
gs = fig.add_gridspec(2, 6)
ax1 = fig.add_subplot(gs[0, 0:2])
ax2 = fig.add_subplot(gs[0, 2:4])
ax3 = fig.add_subplot(gs[0, 4:6])
ax4 = fig.add_subplot(gs[1, 1:3])
ax5 = fig.add_subplot(gs[1, 3:5])
axs = [ax1, ax2, ax3, ax4, ax5]

plt.suptitle('Data yield', fontsize=24)
plt.subplots_adjust(wspace=0.5, hspace=0.5)

for i in range(len(dsdt_index)-1):
    axs[0].errorbar((phi_d_2H_dsdt_minust_high[dsdt_index[i]:dsdt_index[i+1]] + phi_d_2H_dsdt_minust_low[dsdt_index[i]:dsdt_index[i+1]])/2, \
                    phi_d_2H_dsdt_yield_data[0,dsdt_index[i]:dsdt_index[i+1]], \
                    xerr=(phi_d_2H_dsdt_minust_high[dsdt_index[i]:dsdt_index[i+1]] - phi_d_2H_dsdt_minust_low[dsdt_index[i]:dsdt_index[i+1]])/2, \
                    yerr=phi_d_2H_dsdt_yield_data_statserr[0,dsdt_index[i]:dsdt_index[i+1]], \
                    fmt=dsdt_color_list[i]+'.', label=dsdt_label_list[i])
axs[0].set_title(r'$d\sigma / dt$')
axs[0].set_xlabel(r'$-t[GeV^2/c]$')
axs[0].set_ylabel(r'$Y_{\rm data}$')
axs[0].set_xlim(0, 2)
axs[0].set_ylim(0, 400)
axs[0].legend()

# for i in range(len(Wcostheta_index)-1):
#     axs[1].errorbar((phi_d_2H_Wcostheta_costheta_high[Wcostheta_index[i]:Wcostheta_index[i+1]] + phi_d_2H_Wcostheta_costheta_low[Wcostheta_index[i]:Wcostheta_index[i+1]])/2, \
#                     phi_d_2H_Wcostheta_yield_data[0,Wcostheta_index[i]:Wcostheta_index[i+1]], \
#                     xerr=(phi_d_2H_Wcostheta_costheta_high[Wcostheta_index[i]:Wcostheta_index[i+1]] - phi_d_2H_Wcostheta_costheta_low[Wcostheta_index[i]:Wcostheta_index[i+1]])/2, \
#                     yerr=phi_d_2H_Wcostheta_yield_data_statserr[0,Wcostheta_index[i]:Wcostheta_index[i+1]], \
#                     fmt='.', color=W_color_list[i], label=W_label_list[i])
# axs[1].set_title(r'$W(\cos \vartheta)$')
# axs[1].set_xlabel(r'$\cos\vartheta$')
# axs[1].set_ylabel(r'$Y_{\rm data}$')
# axs[1].set_xlim(-1, 1)
# axs[1].set_ylim(0, 1000)
# axs[1].legend()

# for i in range(len(Wdecayphi_index)-1):
#     axs[2].errorbar((phi_d_2H_Wdecayphi_decayphi_high[Wdecayphi_index[i]:Wdecayphi_index[i+1]] + phi_d_2H_Wdecayphi_decayphi_low[Wdecayphi_index[i]:Wdecayphi_index[i+1]])/2, \
#                     phi_d_2H_Wdecayphi_yield_data[0,Wdecayphi_index[i]:Wdecayphi_index[i+1]], \
#                     xerr=(phi_d_2H_Wdecayphi_decayphi_high[Wdecayphi_index[i]:Wdecayphi_index[i+1]] - phi_d_2H_Wdecayphi_decayphi_low[Wdecayphi_index[i]:Wdecayphi_index[i+1]])/2, \
#                     yerr=phi_d_2H_Wdecayphi_yield_data_statserr[0,Wdecayphi_index[i]:Wdecayphi_index[i+1]], \
#                     fmt='.', color=W_color_list[i], label=W_label_list[i])
# axs[2].set_title(r'$W(\varphi)$')
# axs[2].set_xlabel(r'$\varphi$ (deg)')
# axs[2].set_ylabel(r'$Y_{\rm data}$')
# axs[2].set_xlim(-180, 180)
# axs[2].set_ylim(0, 800)
# axs[2].legend()

# for i in range(len(Wpolphi_index)-1):
#     axs[3].errorbar((phi_d_2H_Wpolphi_polphi_high[Wpolphi_index[i]:Wpolphi_index[i+1]] + phi_d_2H_Wpolphi_polphi_low[Wpolphi_index[i]:Wpolphi_index[i+1]])/2, \
#                     phi_d_2H_Wpolphi_yield_data[0,Wpolphi_index[i]:Wpolphi_index[i+1]], \
#                     xerr=(phi_d_2H_Wpolphi_polphi_high[Wpolphi_index[i]:Wpolphi_index[i+1]] - phi_d_2H_Wpolphi_polphi_low[Wpolphi_index[i]:Wpolphi_index[i+1]])/2, \
#                     yerr=phi_d_2H_Wpolphi_yield_data_statserr[0,Wpolphi_index[i]:Wpolphi_index[i+1]], \
#                     fmt='.', color=W_color_list[i], label=W_label_list[i])
# axs[3].set_title(r'$W(\Phi)$')
# axs[3].set_xlabel(r'$\Phi$ (deg)')
# axs[3].set_ylabel(r'$Y_{\rm data}$')
# axs[3].set_xlim(-180, 180)
# axs[3].set_ylim(0, 800)
# axs[3].legend()

# for i in range(len(Wpsi_index)-1):
#     axs[4].errorbar((phi_d_2H_Wpsi_psi_high[Wpsi_index[i]:Wpsi_index[i+1]] + phi_d_2H_Wpsi_psi_low[Wpsi_index[i]:Wpsi_index[i+1]])/2, \
#                     phi_d_2H_Wpsi_yield_data[0,Wpsi_index[i]:Wpsi_index[i+1]], \
#                     xerr=(phi_d_2H_Wpsi_psi_high[Wpsi_index[i]:Wpsi_index[i+1]] - phi_d_2H_Wpsi_psi_low[Wpsi_index[i]:Wpsi_index[i+1]])/2, \
#                     yerr=phi_d_2H_Wpsi_yield_data_statserr[0,Wpsi_index[i]:Wpsi_index[i+1]], \
#                     fmt='.', color=W_color_list[i], label=W_label_list[i])
# axs[4].set_title(r'$W(\psi)$')
# axs[4].set_xlabel(r'$\psi$ (deg)')
# axs[4].set_ylabel(r'$Y_{\rm data}$')
# axs[4].set_xlim(-180, 180)
# axs[4].set_ylim(0, 800)
# axs[4].legend()

file_pdf.savefig()
plt.close()

################################################################# EFFICIENCY #############################################################################################################################

# fig = plt.figure(figsize=(24, 16), dpi=300)
# gs = fig.add_gridspec(2, 6)
# ax1 = fig.add_subplot(gs[0, 0:2])
# ax2 = fig.add_subplot(gs[0, 2:4])
# ax3 = fig.add_subplot(gs[0, 4:6])
# ax4 = fig.add_subplot(gs[1, 1:3])
# ax5 = fig.add_subplot(gs[1, 3:5])
# axs = [ax1, ax2, ax3, ax4, ax5]

# plt.suptitle('Efficiency', fontsize=24)
# plt.subplots_adjust(wspace=0.5, hspace=0.5)

# for i in range(len(dsdt_index)-1):
#     axs[0].errorbar((phi_d_2H_dsdt_minust_high[dsdt_index[i]:dsdt_index[i+1]] + phi_d_2H_dsdt_minust_low[dsdt_index[i]:dsdt_index[i+1]])/2, \
#                     phi_d_2H_dsdt_efficiency[0,dsdt_index[i]:dsdt_index[i+1]], \
#                     xerr=(phi_d_2H_dsdt_minust_high[dsdt_index[i]:dsdt_index[i+1]] - phi_d_2H_dsdt_minust_low[dsdt_index[i]:dsdt_index[i+1]])/2, \
#                     yerr=phi_d_2H_dsdt_efficiency_statserr[0,dsdt_index[i]:dsdt_index[i+1]], \
#                     fmt=dsdt_color_list[i]+'.', label=dsdt_label_list[i])
# axs[0].set_title(r'$d\sigma / dt$')
# axs[0].set_xlabel(r'$-t[GeV^2/c]$')
# axs[0].set_ylabel(r'$\epsilon$')
# axs[0].set_xlim(0, 2)
# axs[0].set_ylim(0, 400)
# axs[0].legend()

# for i in range(len(Wcostheta_index)-1):
#     axs[1].errorbar((phi_d_2H_Wcostheta_costheta_high[Wcostheta_index[i]:Wcostheta_index[i+1]] + phi_d_2H_Wcostheta_costheta_low[Wcostheta_index[i]:Wcostheta_index[i+1]])/2, \
#                     phi_d_2H_Wcostheta_efficiency[0,Wcostheta_index[i]:Wcostheta_index[i+1]], \
#                     xerr=(phi_d_2H_Wcostheta_costheta_high[Wcostheta_index[i]:Wcostheta_index[i+1]] - phi_d_2H_Wcostheta_costheta_low[Wcostheta_index[i]:Wcostheta_index[i+1]])/2, \
#                     yerr=phi_d_2H_Wcostheta_efficiency_statserr[0,Wcostheta_index[i]:Wcostheta_index[i+1]], \
#                     fmt='.', color=W_color_list[i], label=W_label_list[i])
# axs[1].set_title(r'$W(\cos \vartheta)$')
# axs[1].set_xlabel(r'$\cos\vartheta$')
# axs[1].set_ylabel(r'$\epsilon$')
# axs[1].set_xlim(-1, 1)
# axs[1].set_ylim(0, 1000)
# axs[1].legend()

# for i in range(len(Wdecayphi_index)-1):
#     axs[2].errorbar((phi_d_2H_Wdecayphi_decayphi_high[Wdecayphi_index[i]:Wdecayphi_index[i+1]] + phi_d_2H_Wdecayphi_decayphi_low[Wdecayphi_index[i]:Wdecayphi_index[i+1]])/2, \
#                     phi_d_2H_Wdecayphi_efficiency[0,Wdecayphi_index[i]:Wdecayphi_index[i+1]], \
#                     xerr=(phi_d_2H_Wdecayphi_decayphi_high[Wdecayphi_index[i]:Wdecayphi_index[i+1]] - phi_d_2H_Wdecayphi_decayphi_low[Wdecayphi_index[i]:Wdecayphi_index[i+1]])/2, \
#                     yerr=phi_d_2H_Wdecayphi_efficiency_statserr[0,Wdecayphi_index[i]:Wdecayphi_index[i+1]], \
#                     fmt='.', color=W_color_list[i], label=W_label_list[i])
# axs[2].set_title(r'$W(\varphi)$')
# axs[2].set_xlabel(r'$\varphi$ (deg)')
# axs[2].set_ylabel(r'$\epsilon$')
# axs[2].set_xlim(-180, 180)
# axs[2].set_ylim(0, 800)
# axs[2].legend()

# for i in range(len(Wpolphi_index)-1):
#     axs[3].errorbar((phi_d_2H_Wpolphi_polphi_high[Wpolphi_index[i]:Wpolphi_index[i+1]] + phi_d_2H_Wpolphi_polphi_low[Wpolphi_index[i]:Wpolphi_index[i+1]])/2, \
#                     phi_d_2H_Wpolphi_efficiency[0,Wpolphi_index[i]:Wpolphi_index[i+1]], \
#                     xerr=(phi_d_2H_Wpolphi_polphi_high[Wpolphi_index[i]:Wpolphi_index[i+1]] - phi_d_2H_Wpolphi_polphi_low[Wpolphi_index[i]:Wpolphi_index[i+1]])/2, \
#                     yerr=phi_d_2H_Wpolphi_efficiency_statserr[0,Wpolphi_index[i]:Wpolphi_index[i+1]], \
#                     fmt='.', color=W_color_list[i], label=W_label_list[i])
# axs[3].set_title(r'$W(\Phi)$')
# axs[3].set_xlabel(r'$\Phi$ (deg)')
# axs[3].set_ylabel(r'$\epsilon$')
# axs[3].set_xlim(-180, 180)
# axs[3].set_ylim(0, 800)
# axs[3].legend()

# for i in range(len(Wpsi_index)-1):
#     axs[4].errorbar((phi_d_2H_Wpsi_psi_high[Wpsi_index[i]:Wpsi_index[i+1]] + phi_d_2H_Wpsi_psi_low[Wpsi_index[i]:Wpsi_index[i+1]])/2, \
#                     phi_d_2H_Wpsi_efficiency[0,Wpsi_index[i]:Wpsi_index[i+1]], \
#                     xerr=(phi_d_2H_Wpsi_psi_high[Wpsi_index[i]:Wpsi_index[i+1]] - phi_d_2H_Wpsi_psi_low[Wpsi_index[i]:Wpsi_index[i+1]])/2, \
#                     yerr=phi_d_2H_Wpsi_efficiency_statserr[0,Wpsi_index[i]:Wpsi_index[i+1]], \
#                     fmt='.', color=W_color_list[i], label=W_label_list[i])
# axs[4].set_title(r'$W(\psi)$')
# axs[4].set_xlabel(r'$\psi$ (deg)')
# axs[4].set_ylabel(r'$\epsilon$')
# axs[4].set_xlim(-180, 180)
# axs[4].set_ylim(0, 800)
# axs[4].legend()

# file_pdf.savefig()
# plt.close()

################################################################# SYSTEMATICS: DATA YIELD #################################################################################################################

for k in range(len(syst_list)):
    if(syst_list[k] == 'Simulation weight'):
        continue

    fig = plt.figure(figsize=(24, 16), dpi=300)
    gs = fig.add_gridspec(2, 6)
    ax1 = fig.add_subplot(gs[0, 0:2])
    ax2 = fig.add_subplot(gs[0, 2:4])
    ax3 = fig.add_subplot(gs[0, 4:6])
    ax4 = fig.add_subplot(gs[1, 1:3])
    ax5 = fig.add_subplot(gs[1, 3:5])
    axs = [ax1, ax2, ax3, ax4, ax5]

    plt.suptitle('Data yield relative difference: '+syst_list[k],fontsize=24)
    plt.subplots_adjust(wspace=0.5, hspace=0.5)

    for i in range(1, phi_d_2H_dsdt_yield_data_statserr.shape[0]):
        if (subplot_list[i] == k):
            for j in range(len(dsdt_index)-1):
                axs[0].scatter( phi_d_2H_dsdt_minust_center[dsdt_index[j]:dsdt_index[j+1]], \
                                (phi_d_2H_dsdt_yield_data[0,dsdt_index[j]:dsdt_index[j+1]]-phi_d_2H_dsdt_yield_data[i,dsdt_index[j]:dsdt_index[j+1]])/(phi_d_2H_dsdt_yield_data[0,dsdt_index[j]:dsdt_index[j+1]]), \
                                label=tag_list[i], color=color_list[i])

    # for i in range(1, phi_d_2H_Wcostheta_yield_data_statserr.shape[0]):
    #     if (subplot_list[i] == k):
    #         for j in range(len(Wcostheta_index)-1):
    #             axs[1].scatter( phi_d_2H_Wcostheta_costheta_center[Wcostheta_index[j]:Wcostheta_index[j+1]], \
    #                             (phi_d_2H_Wcostheta_yield_data[0,Wcostheta_index[j]:Wcostheta_index[j+1]]-phi_d_2H_Wcostheta_yield_data[i,Wcostheta_index[j]:Wcostheta_index[j+1]])/(phi_d_2H_Wcostheta_yield_data[0,Wcostheta_index[j]:Wcostheta_index[j+1]]), \
    #                             label=tag_list[i], color=color_list[i])

    # for i in range(1, phi_d_2H_Wdecayphi_yield_data_statserr.shape[0]):
    #     if (subplot_list[i] == k):
    #         for j in range(len(Wdecayphi_index)-1):
    #             axs[2].scatter( phi_d_2H_Wdecayphi_decayphi_center[Wdecayphi_index[j]:Wdecayphi_index[j+1]], \
    #                             (phi_d_2H_Wdecayphi_yield_data[0,Wdecayphi_index[j]:Wdecayphi_index[j+1]]-phi_d_2H_Wdecayphi_yield_data[i,Wdecayphi_index[j]:Wdecayphi_index[j+1]])/(phi_d_2H_Wdecayphi_yield_data[0,Wdecayphi_index[j]:Wdecayphi_index[j+1]]), \
    #                             label=tag_list[i], color=color_list[i])

    # for i in range(1, phi_d_2H_Wpolphi_yield_data_statserr.shape[0]):
    #     if (subplot_list[i] == k):
    #         for j in range(len(Wpolphi_index)-1):
    #             axs[3].scatter( phi_d_2H_Wpolphi_polphi_center[Wpolphi_index[j]:Wpolphi_index[j+1]], \
    #                             (phi_d_2H_Wpolphi_yield_data[0,Wpolphi_index[j]:Wpolphi_index[j+1]]-phi_d_2H_Wpolphi_yield_data[i,Wpolphi_index[j]:Wpolphi_index[j+1]])/(phi_d_2H_Wpolphi_yield_data[0,Wpolphi_index[j]:Wpolphi_index[j+1]]), \
    #                             label=tag_list[i], color=color_list[i])

    # for i in range(1, phi_d_2H_Wpsi_yield_data_statserr.shape[0]):
    #     if (subplot_list[i] == k):
    #         for j in range(len(Wpsi_index)-1):
    #             axs[4].scatter( phi_d_2H_Wpsi_psi_center[Wpsi_index[j]:Wpsi_index[j+1]], \
    #                             (phi_d_2H_Wpsi_yield_data[0,Wpsi_index[j]:Wpsi_index[j+1]]-phi_d_2H_Wpsi_yield_data[i,Wpsi_index[j]:Wpsi_index[j+1]])/(phi_d_2H_Wpsi_yield_data[0,Wpsi_index[j]:Wpsi_index[j+1]]), \
    #                             label=tag_list[i], color=color_list[i])

    for i, ax in enumerate(axs):
        ax.set_xlim([[0, -1, -180, -180, -180][i], [2, 1, 180, 180, 180][i]])
        ax.set_ylim(-0.5, 0.5)
        ax.set_xlabel(xlabel_list[i])
        ax.set_ylabel(r'$(Y^\prime_{\rm data} - Y_{\rm data})/Y_{\rm data}$')
        ax.set_title(title_list[i])
        ax.legend()
        legend_without_duplicate_labels(ax)
        ax.plot([ax.get_xlim()[0], ax.get_xlim()[1]], [0.1, 0.1], 'r--')
        ax.plot([ax.get_xlim()[0], ax.get_xlim()[1]], [-0.1, -0.1], 'r--')

    file_pdf.savefig()
    plt.close()

################################################################# SYSTEMATICS: SIM YIELD ####################################################################################################################

# for k in range(len(syst_list)):
#     if(syst_list[k] == 'Simulation weight'):
#         continue

#     fig = plt.figure(figsize=(24, 16), dpi=300)
#     gs = fig.add_gridspec(2, 6)
#     ax1 = fig.add_subplot(gs[0, 0:2])
#     ax2 = fig.add_subplot(gs[0, 2:4])
#     ax3 = fig.add_subplot(gs[0, 4:6])
#     ax4 = fig.add_subplot(gs[1, 1:3])
#     ax5 = fig.add_subplot(gs[1, 3:5])
#     axs = [ax1, ax2, ax3, ax4, ax5]

#     plt.suptitle('Simulation yield relative difference: '+syst_list[k],fontsize=24)
#     plt.subplots_adjust(wspace=0.5, hspace=0.5)

#     for i in range(1, phi_d_2H_dsdt_yield_sim_statserr.shape[0]):
#         print(i)
#         if (subplot_list[i] == k):
#             for j in range(len(dsdt_index)-1):
#                 axs[0].scatter( phi_d_2H_dsdt_minust_center[dsdt_index[j]:dsdt_index[j+1]], \
#                                 (phi_d_2H_dsdt_yield_sim[0,dsdt_index[j]:dsdt_index[j+1]]-phi_d_2H_dsdt_yield_sim[i,dsdt_index[j]:dsdt_index[j+1]])/(phi_d_2H_dsdt_yield_sim[0,dsdt_index[j]:dsdt_index[j+1]]), \
#                                 label=tag_list[i], color=color_list[i])

#     for i in range(1, phi_d_2H_Wcostheta_yield_sim_statserr.shape[0]):
#         if (subplot_list[i] == k):
#             for j in range(len(Wcostheta_index)-1):
#                 axs[1].scatter( phi_d_2H_Wcostheta_costheta_center[Wcostheta_index[j]:Wcostheta_index[j+1]], \
#                                 (phi_d_2H_Wcostheta_yield_sim[0,Wcostheta_index[j]:Wcostheta_index[j+1]]-phi_d_2H_Wcostheta_yield_sim[i,Wcostheta_index[j]:Wcostheta_index[j+1]])/(phi_d_2H_Wcostheta_yield_sim[0,Wcostheta_index[j]:Wcostheta_index[j+1]]), \
#                                 label=tag_list[i], color=color_list[i])

#     for i in range(1, phi_d_2H_Wdecayphi_yield_sim_statserr.shape[0]):
#         if (subplot_list[i] == k):
#             for j in range(len(Wdecayphi_index)-1):
#                 axs[2].scatter( phi_d_2H_Wdecayphi_decayphi_center[Wdecayphi_index[j]:Wdecayphi_index[j+1]], \
#                                 (phi_d_2H_Wdecayphi_yield_sim[0,Wdecayphi_index[j]:Wdecayphi_index[j+1]]-phi_d_2H_Wdecayphi_yield_sim[i,Wdecayphi_index[j]:Wdecayphi_index[j+1]])/(phi_d_2H_Wdecayphi_yield_sim[0,Wdecayphi_index[j]:Wdecayphi_index[j+1]]), \
#                                 label=tag_list[i], color=color_list[i])

#     for i in range(1, phi_d_2H_Wpolphi_yield_sim_statserr.shape[0]):
#         if (subplot_list[i] == k):
#             for j in range(len(Wpolphi_index)-1):
#                 axs[3].scatter( phi_d_2H_Wpolphi_polphi_center[Wpolphi_index[j]:Wpolphi_index[j+1]], \
#                                 (phi_d_2H_Wpolphi_yield_sim[0,Wpolphi_index[j]:Wpolphi_index[j+1]]-phi_d_2H_Wpolphi_yield_sim[i,Wpolphi_index[j]:Wpolphi_index[j+1]])/(phi_d_2H_Wpolphi_yield_sim[0,Wpolphi_index[j]:Wpolphi_index[j+1]]), \
#                                 label=tag_list[i], color=color_list[i])

#     for i in range(1, phi_d_2H_Wpsi_yield_sim_statserr.shape[0]):
#         if (subplot_list[i] == k):
#             for j in range(len(Wpsi_index)-1):
#                 axs[4].scatter( phi_d_2H_Wpsi_psi_center[Wpsi_index[j]:Wpsi_index[j+1]], \
#                                 (phi_d_2H_Wpsi_yield_sim[0,Wpsi_index[j]:Wpsi_index[j+1]]-phi_d_2H_Wpsi_yield_sim[i,Wpsi_index[j]:Wpsi_index[j+1]])/(phi_d_2H_Wpsi_yield_sim[0,Wpsi_index[j]:Wpsi_index[j+1]]), \
#                                 label=tag_list[i], color=color_list[i])

#     for i, ax in enumerate(axs):
#         ax.set_xlim([[0, -1, -180, -180, -180][i], [2, 1, 180, 180, 180][i]])
#         ax.set_ylim(-1, 1)
#         ax.set_xlabel(xlabel_list[i])
#         ax.set_ylabel(r'$(Y^\prime_{\rm sim} - Y_{\rm sim})/Y_{\rm sim}$')
#         ax.set_title(title_list[i])
#         ax.legend()
#         legend_without_duplicate_labels(ax)
#         ax.plot([ax.get_xlim()[0], ax.get_xlim()[1]], [0.1, 0.1], 'r--')
#         ax.plot([ax.get_xlim()[0], ax.get_xlim()[1]], [-0.1, -0.1], 'r--')

#     file_pdf.savefig()
#     plt.close()

################################################################# SYSTEMATICS: BARLOW SCORE ###############################################################################################################

# for k in range(len(syst_list)):
#     if(syst_list[k] == 'Simulation weight'):
#         continue

#     fig = plt.figure(figsize=(24, 16), dpi=300)
#     gs = fig.add_gridspec(2, 6)
#     ax1 = fig.add_subplot(gs[0, 0:2])
#     ax2 = fig.add_subplot(gs[0, 2:4])
#     ax3 = fig.add_subplot(gs[0, 4:6])
#     ax4 = fig.add_subplot(gs[1, 1:3])
#     ax5 = fig.add_subplot(gs[1, 3:5])
#     axs = [ax1, ax2, ax3, ax4, ax5]

#     plt.suptitle('Barlow score: '+syst_list[k],fontsize=24)
#     plt.subplots_adjust(wspace=0.5, hspace=0.5)

#     for i in range(1, phi_d_2H_dsdt_yield_sim_statserr.shape[0]):
#         if (subplot_list[i] == k):
#             for j in range(len(dsdt_index)-1):
#                axs[0].scatter( phi_d_2H_dsdt_minust_center[dsdt_index[j]:dsdt_index[j+1]], \
#                                (phi_d_2H_dsdt_results[0,dsdt_index[j]:dsdt_index[j+1]] - phi_d_2H_dsdt_results[i,dsdt_index[j]:dsdt_index[j+1]])/(1E-40+np.sqrt(np.abs(phi_d_2H_dsdt_results_statserr[0,dsdt_index[j]:dsdt_index[j+1]]**2 - phi_d_2H_dsdt_results_statserr[i,dsdt_index[j]:dsdt_index[j+1]]**2))), \
#                                label=tag_list[i], color=color_list[i])

#     for i in range(1, phi_d_2H_Wcostheta_yield_sim_statserr.shape[0]):
#         if (subplot_list[i] == k):
#             for j in range(len(Wcostheta_index)-1):
#                 axs[1].scatter( phi_d_2H_Wcostheta_costheta_center[Wcostheta_index[j]:Wcostheta_index[j+1]], \
#                                 (phi_d_2H_Wcostheta_results[0,Wcostheta_index[j]:Wcostheta_index[j+1]] - phi_d_2H_Wcostheta_results[i,Wcostheta_index[j]:Wcostheta_index[j+1]])/(1E-40+np.sqrt(np.abs(phi_d_2H_Wcostheta_results_statserr[0,Wcostheta_index[j]:Wcostheta_index[j+1]]**2 - phi_d_2H_Wcostheta_results_statserr[i,Wcostheta_index[j]:Wcostheta_index[j+1]]**2))), \
#                                 label=tag_list[i], color=color_list[i])

#     for i in range(1, phi_d_2H_Wdecayphi_yield_sim_statserr.shape[0]):
#         if (subplot_list[i] == k):
#             for j in range(len(Wdecayphi_index)-1):
#                 axs[2].scatter( phi_d_2H_Wdecayphi_decayphi_center[Wdecayphi_index[j]:Wdecayphi_index[j+1]], \
#                                 (phi_d_2H_Wdecayphi_results[0,Wdecayphi_index[j]:Wdecayphi_index[j+1]] - phi_d_2H_Wdecayphi_results[i,Wdecayphi_index[j]:Wdecayphi_index[j+1]])/(1E-40+np.sqrt(np.abs(phi_d_2H_Wdecayphi_results_statserr[0,Wdecayphi_index[j]:Wdecayphi_index[j+1]]**2 - phi_d_2H_Wdecayphi_results_statserr[i,Wdecayphi_index[j]:Wdecayphi_index[j+1]]**2))), \
#                                 label=tag_list[i], color=color_list[i])

#     for i in range(1, phi_d_2H_Wpolphi_yield_sim_statserr.shape[0]):
#         if (subplot_list[i] == k):
#             for j in range(len(Wpolphi_index)-1):
#                 axs[3].scatter( phi_d_2H_Wpolphi_polphi_center[Wpolphi_index[j]:Wpolphi_index[j+1]], \
#                                 (phi_d_2H_Wpolphi_results[0,Wpolphi_index[j]:Wpolphi_index[j+1]] - phi_d_2H_Wpolphi_results[i,Wpolphi_index[j]:Wpolphi_index[j+1]])/(1E-40+np.sqrt(np.abs(phi_d_2H_Wpolphi_results_statserr[0,Wpolphi_index[j]:Wpolphi_index[j+1]]**2 - phi_d_2H_Wpolphi_results_statserr[i,Wpolphi_index[j]:Wpolphi_index[j+1]]**2))), \
#                                 label=tag_list[i], color=color_list[i])

#     for i in range(1, phi_d_2H_Wpsi_yield_sim_statserr.shape[0]):
#         if (subplot_list[i] == k):
#             for j in range(len(Wpsi_index)-1):
#                 axs[4].scatter( phi_d_2H_Wpsi_psi_center[Wpsi_index[j]:Wpsi_index[j+1]], \
#                                 (phi_d_2H_Wpsi_results[0,Wpsi_index[j]:Wpsi_index[j+1]]-phi_d_2H_Wpsi_results[i,Wpsi_index[j]:Wpsi_index[j+1]])/(1E-40+np.sqrt(np.abs(phi_d_2H_Wpsi_results_statserr[0,Wpsi_index[j]:Wpsi_index[j+1]]**2 - phi_d_2H_Wpsi_results_statserr[i,Wpsi_index[j]:Wpsi_index[j+1]]**2))), \
#                                 label=tag_list[i], color=color_list[i])

#     for i, ax in enumerate(axs):
#         ax.set_xlim([[0, -1, -180, -180, -180][i], [2, 1, 180, 180, 180][i]])
#         ax.set_ylim(-15, 15)
#         ax.set_xlabel(xlabel_list[i])
#         ax.set_ylabel(r'Barlow score')
#         ax.set_title(title_list[i])
#         ax.legend()
#         legend_without_duplicate_labels(ax)
#         ax.fill_between([0, 2], [1, 1], [-1, -1], color='green', alpha=0.1)
#         ax.fill_between([0, 2], [4, 4], [1, 1], color='yellow', alpha=0.1)
#         ax.fill_between([0, 2], [15, 15], [4, 4], color='red', alpha=0.1)
#         ax.fill_between([0, 2], [-1, -1], [-4, -4], color='yellow', alpha=0.1)
#         ax.fill_between([0, 2], [-4, -4], [-15, -15], color='red', alpha=0.1)

#     file_pdf.savefig()
#     plt.close()

################################################################# SYSTEMATICS: OBSERVABLE #################################################################################################################

# for k in range(len(syst_list)):
#     if(syst_list[k] == 'Simulation weight'):
#         continue

#     fig = plt.figure(figsize=(24, 16), dpi=300)
#     gs = fig.add_gridspec(2, 6)
#     ax1 = fig.add_subplot(gs[0, 0:2])
#     ax2 = fig.add_subplot(gs[0, 2:4])
#     ax3 = fig.add_subplot(gs[0, 4:6])
#     ax4 = fig.add_subplot(gs[1, 1:3])
#     ax5 = fig.add_subplot(gs[1, 3:5])
#     axs = [ax1, ax2, ax3, ax4, ax5]

#     plt.suptitle('Observable relative difference: '+syst_list[k],fontsize=24)
#     plt.subplots_adjust(wspace=0.5, hspace=0.5)

#     for i in range(1, phi_d_2H_dsdt_yield_sim_statserr.shape[0]):
#         if (subplot_list[i] == k):
#             for j in range(len(dsdt_index)-1):
#                 axs[0].scatter( phi_d_2H_dsdt_minust_center[dsdt_index[j]:dsdt_index[j+1]], \
#                                 (phi_d_2H_dsdt_results[0,dsdt_index[j]:dsdt_index[j+1]] - phi_d_2H_dsdt_results[i,dsdt_index[j]:dsdt_index[j+1]])/phi_d_2H_dsdt_results[0,dsdt_index[j]:dsdt_index[j+1]], \
#                                 label=tag_list[i], color=color_list[i])

#     for i in range(1, phi_d_2H_Wcostheta_yield_sim_statserr.shape[0]):
#         if (subplot_list[i] == k):
#             for j in range(len(Wcostheta_index)-1):
#                 axs[1].scatter( phi_d_2H_Wcostheta_costheta_center[Wcostheta_index[j]:Wcostheta_index[j+1]], \
#                                 (phi_d_2H_Wcostheta_results[0,Wcostheta_index[j]:Wcostheta_index[j+1]] - phi_d_2H_Wcostheta_results[i,Wcostheta_index[j]:Wcostheta_index[j+1]])/phi_d_2H_Wcostheta_results[0,Wcostheta_index[j]:Wcostheta_index[j+1]], \
#                                 label=tag_list[i], color=color_list[i])

#     for i in range(1, phi_d_2H_Wdecayphi_yield_sim_statserr.shape[0]):
#         if (subplot_list[i] == k):
#             for j in range(len(Wdecayphi_index)-1):
#                 axs[2].scatter( phi_d_2H_Wdecayphi_decayphi_center[Wdecayphi_index[j]:Wdecayphi_index[j+1]], \
#                                 (phi_d_2H_Wdecayphi_results[0,Wdecayphi_index[j]:Wdecayphi_index[j+1]] - phi_d_2H_Wdecayphi_results[i,Wdecayphi_index[j]:Wdecayphi_index[j+1]])/phi_d_2H_Wdecayphi_results[0,Wdecayphi_index[j]:Wdecayphi_index[j+1]], \
#                                 label=tag_list[i], color=color_list[i])

#     for i in range(1, phi_d_2H_Wpolphi_yield_sim_statserr.shape[0]):
#         if (subplot_list[i] == k):
#             for j in range(len(Wpolphi_index)-1):
#                 axs[3].scatter( phi_d_2H_Wpolphi_polphi_center[Wpolphi_index[j]:Wpolphi_index[j+1]], \
#                                 (phi_d_2H_Wpolphi_results[0,Wpolphi_index[j]:Wpolphi_index[j+1]] - phi_d_2H_Wpolphi_results[i,Wpolphi_index[j]:Wpolphi_index[j+1]])/phi_d_2H_Wpolphi_results[0,Wpolphi_index[j]:Wpolphi_index[j+1]], \
#                                 label=tag_list[i], color=color_list[i])

#     for i in range(1, phi_d_2H_Wpsi_yield_sim_statserr.shape[0]):
#         if (subplot_list[i] == k):
#             for j in range(len(Wpsi_index)-1):
#                 axs[4].scatter( phi_d_2H_Wpsi_psi_center[Wpsi_index[j]:Wpsi_index[j+1]], \
#                                 (phi_d_2H_Wpsi_results[0,Wpsi_index[j]:Wpsi_index[j+1]]-phi_d_2H_Wpsi_results[i,Wpsi_index[j]:Wpsi_index[j+1]])/phi_d_2H_Wpsi_results[0,Wpsi_index[j]:Wpsi_index[j+1]], \
#                                 label=tag_list[i], color=color_list[i])

#     for i, ax in enumerate(axs):
#         ax.set_xlim([[0, -1, -180, -180, -180][i], [2, 1, 180, 180, 180][i]])
#         ax.set_ylim(-0.2, 0.2)
#         ax.set_xlabel(xlabel_list[i])
#         ax.set_ylabel(r'$(\sigma^\prime - \sigma)/\sigma$')
#         ax.set_title(title_list[i])
#         ax.legend()
#         legend_without_duplicate_labels(ax)
#         ax.fill_between([0, 2], [1, 1], [-1, -1], color='green', alpha=0.1)
#         ax.fill_between([0, 2], [4, 4], [1, 1], color='yellow', alpha=0.1)
#         ax.fill_between([0, 2], [15, 15], [4, 4], color='red', alpha=0.1)
#         ax.fill_between([0, 2], [-1, -1], [-4, -4], color='yellow', alpha=0.1)
#         ax.fill_between([0, 2], [-4, -4], [-15, -15], color='red', alpha=0.1)

#     file_pdf.savefig()
#     plt.close()

################################################################# SYSTEMATICS: UNCERTAINTIES #############################################################################################################

# fig = plt.figure(figsize=(6*num_col, 6*(num_row+1)), dpi=300)
# gs = fig.add_gridspec(num_row+1, num_col)
# axs = gs.subplots()
# total_uncertainties = [np.zeros_like(phi_d_2H_dsdt_results[0,dsdt_index[j]:dsdt_index[j+1]]) for j in range(len(dsdt_index)-1)]
# for j in range(len(dsdt_index)-1):
#     for i in range(1, phi_d_2H_dsdt_yield_data_statserr.shape[0]):
#         if i == 1:
#             temp_array = (phi_d_2H_dsdt_results[0,dsdt_index[j]:dsdt_index[j+1]] - phi_d_2H_dsdt_results[2,dsdt_index[j]:dsdt_index[j+1]])/phi_d_2H_dsdt_results[0,dsdt_index[j]:dsdt_index[j+1]]
#         elif subplot_list[i] != subplot_list[i-1]:
#             subplot_row = subplot_list[i-1]//num_col
#             subplot_col = subplot_list[i-1]%num_col
#             temp_std = np.std(temp_array, axis=0)
#             axs[subplot_row, subplot_col].scatter(phi_d_2H_dsdt_minust_center[dsdt_index[j]:dsdt_index[j+1]], temp_std, label=legend_list[j])
#             total_uncertainties[j] += temp_std**2
#             temp_array = (phi_d_2H_dsdt_results[0,dsdt_index[j]:dsdt_index[j+1]] - phi_d_2H_dsdt_results[2,dsdt_index[j]:dsdt_index[j+1]])/phi_d_2H_dsdt_results[0,dsdt_index[j]:dsdt_index[j+1]]
#         if i == 1:
#             continue
#         temp_array = np.vstack((temp_array, (phi_d_2H_dsdt_results[0,dsdt_index[j]:dsdt_index[j+1]] - phi_d_2H_dsdt_results[i,dsdt_index[j]:dsdt_index[j+1]])/phi_d_2H_dsdt_results[0,dsdt_index[j]:dsdt_index[j+1]]))
#         if i == phi_d_2H_dsdt_yield_data_statserr.shape[0]-1:
#             subplot_row = subplot_list[i-1]//num_col
#             subplot_col = subplot_list[i-1]%num_col
#             temp_std = np.std(temp_array, axis=0)
#             axs[subplot_row, subplot_col].scatter(phi_d_2H_dsdt_minust_center[dsdt_index[j]:dsdt_index[j+1]], temp_std, label=legend_list[j])
#             total_uncertainties[j] += temp_std**2
# for j in range(len(dsdt_index)-1):
#     total_uncertainties[j] = np.sqrt(total_uncertainties[j])
#     axs[4, 0].scatter(phi_d_2H_dsdt_minust_center[dsdt_index[j]:dsdt_index[j+1]], total_uncertainties[j])
#     axs[4, 0].legend()
#     axs[4, 0].set_xlabel(r'$-t[GeV^2/c]$')
#     axs[4, 0].set_ylabel(r'$\delta\sigma/\sigma$')
#     axs[4, 0].set_xlim(0, 2)
#     axs[4, 0].set_ylim(0, 0.2)
#     axs[4, 0].set_title('Total point-to-point uncertainty')

# for subplot_row in range(num_row):
#     for subplot_col in range(num_col):
#         axs[subplot_row, subplot_col].legend()
#         axs[subplot_row, subplot_col].set_xlabel(r'$-t[GeV^2/c]$')
#         axs[subplot_row, subplot_col].set_ylabel(r'$\delta\sigma/\sigma$')
#         axs[subplot_row, subplot_col].set_xlim(0, 2)
#         axs[subplot_row, subplot_col].set_ylim(0, 0.2)
#         axs[subplot_row, subplot_col].set_title(title_list[subplot_list[i]])
# plt.suptitle('Observable uncertainties')
# file_pdf.savefig()
# plt.close()




# # for i in range(8):
# #     subplot_row = i//num_col
# #     subplot_col = i%num_col
# #     for j in range(len(dsdt_index)-1):
# #         if (i == 0):
# #             temp_array = (phi_d_2H_dsdt_results[0,dsdt_index[j]:dsdt_index[j+1]] - phi_d_2H_dsdt_results[2,dsdt_index[j]:dsdt_index[j+1]])/phi_d_2H_dsdt_results[0,dsdt_index[j]:dsdt_index[j+1]]
# #             temp_array = np.vstack((temp_array, (phi_d_2H_dsdt_results[0,dsdt_index[j]:dsdt_index[j+1]] - phi_d_2H_dsdt_results[3,dsdt_index[j]:dsdt_index[j+1]])/phi_d_2H_dsdt_results[0,dsdt_index[j]:dsdt_index[j+1]]))
# #             temp_array = np.vstack((temp_array, (phi_d_2H_dsdt_results[0,dsdt_index[j]:dsdt_index[j+1]] - phi_d_2H_dsdt_results[4,dsdt_index[j]:dsdt_index[j+1]])/phi_d_2H_dsdt_results[0,dsdt_index[j]:dsdt_index[j+1]]))
# #             temp_array = np.vstack((temp_array, np.zeros_like(phi_d_2H_dsdt_results[0,dsdt_index[j]:dsdt_index[j+1]])))
# #         elif (i == 1):
# #             temp_array = (phi_d_2H_dsdt_results[0,dsdt_index[j]:dsdt_index[j+1]] - phi_d_2H_dsdt_results[5,dsdt_index[j]:dsdt_index[j+1]])/phi_d_2H_dsdt_results[0,dsdt_index[j]:dsdt_index[j+1]]
# #             for k in range(1,25):
# #                 temp_array = np.vstack((temp_array, (phi_d_2H_dsdt_results[0,dsdt_index[j]:dsdt_index[j+1]] - phi_d_2H_dsdt_results[5+k,dsdt_index[j]:dsdt_index[j+1]])/phi_d_2H_dsdt_results[0,dsdt_index[j]:dsdt_index[j+1]]))
# #         elif (i == 6):
# #             temp_array = (phi_d_2H_dsdt_results[0,dsdt_index[j]:dsdt_index[j+1]] - phi_d_2H_dsdt_results[46,dsdt_index[j]:dsdt_index[j+1]])/phi_d_2H_dsdt_results[0,dsdt_index[j]:dsdt_index[j+1]]
# #             temp_array = np.vstack((temp_array, (phi_d_2H_dsdt_results[0,dsdt_index[j]:dsdt_index[j+1]] - phi_d_2H_dsdt_results[47,dsdt_index[j]:dsdt_index[j+1]])/phi_d_2H_dsdt_results[0,dsdt_index[j]:dsdt_index[j+1]]))
# #             temp_array = np.vstack((temp_array, (phi_d_2H_dsdt_results[0,dsdt_index[j]:dsdt_index[j+1]] - phi_d_2H_dsdt_results[48,dsdt_index[j]:dsdt_index[j+1]])/phi_d_2H_dsdt_results[0,dsdt_index[j]:dsdt_index[j+1]]))
# #             temp_array = np.vstack((temp_array, np.zeros_like(phi_d_2H_dsdt_results[0,dsdt_index[j]:dsdt_index[j+1]])))
# #         elif (i == 7):
# #             temp_array = (phi_d_2H_dsdt_results[0,dsdt_index[j]:dsdt_index[j+1]] - phi_d_2H_dsdt_results[49,dsdt_index[j]:dsdt_index[j+1]])/phi_d_2H_dsdt_results[0,dsdt_index[j]:dsdt_index[j+1]]
# #             # temp_array = np.vstack((temp_array, (phi_d_2H_dsdt_results[0,dsdt_index[j]:dsdt_index[j+1]] - phi_d_2H_dsdt_results[50,dsdt_index[j]:dsdt_index[j+1]])/phi_d_2H_dsdt_results[0,dsdt_index[j]:dsdt_index[j+1]]))
# #             temp_array = np.vstack((temp_array, np.zeros_like(phi_d_2H_dsdt_results[0,dsdt_index[j]:dsdt_index[j+1]])))
# #         else:
# #             temp_array = (phi_d_2H_dsdt_results[0,dsdt_index[j]:dsdt_index[j+1]] - phi_d_2H_dsdt_results[30+(i-2)*4+0,dsdt_index[j]:dsdt_index[j+1]])/phi_d_2H_dsdt_results[0,dsdt_index[j]:dsdt_index[j+1]]
# #             temp_array = np.vstack((temp_array, (phi_d_2H_dsdt_results[0,dsdt_index[j]:dsdt_index[j+1]] - phi_d_2H_dsdt_results[30+(i-2)*4+1,dsdt_index[j]:dsdt_index[j+1]])/phi_d_2H_dsdt_results[0,dsdt_index[j]:dsdt_index[j+1]]))
# #             temp_array = np.vstack((temp_array, (phi_d_2H_dsdt_results[0,dsdt_index[j]:dsdt_index[j+1]] - phi_d_2H_dsdt_results[30+(i-2)*4+2,dsdt_index[j]:dsdt_index[j+1]])/phi_d_2H_dsdt_results[0,dsdt_index[j]:dsdt_index[j+1]]))
# #             temp_array = np.vstack((temp_array, (phi_d_2H_dsdt_results[0,dsdt_index[j]:dsdt_index[j+1]] - phi_d_2H_dsdt_results[30+(i-2)*4+3,dsdt_index[j]:dsdt_index[j+1]])/phi_d_2H_dsdt_results[0,dsdt_index[j]:dsdt_index[j+1]]))
# #             temp_array = np.vstack((temp_array, np.zeros_like(phi_d_2H_dsdt_results[0,dsdt_index[j]:dsdt_index[j+1]])))
# #         temp_std = np.std(temp_array, axis=0)
# #         axs[subplot_row, subplot_col].scatter(phi_d_2H_dsdt_minust_center[dsdt_index[j]:dsdt_index[j+1]], temp_std, label=legend_list[j])
# #         total_uncertainties[j] += temp_std**2
# #     # axs[subplot_row, subplot_col].legend()
# #     axs[subplot_row, subplot_col].set_xlabel(r'$-t[GeV^2/c]$')
# #     axs[subplot_row, subplot_col].set_ylabel(r'$\delta\sigma/\sigma$')
# #     axs[subplot_row, subplot_col].set_xlim(0, 2)
# #     axs[subplot_row, subplot_col].set_ylim(0, 0.2)
# #     axs[subplot_row, subplot_col].set_title(title_list[i])
# # for j in range(len(index)-1):
# #     total_uncertainties[j] = np.sqrt(total_uncertainties[j])
# #     axs[2, 2].scatter(phi_d_2H_dsdt_minust_center[dsdt_index[j]:dsdt_index[j+1]], total_uncertainties[j], label=legend_list[j])
# #     axs[2, 2].legend()
# #     axs[2, 2].set_xlabel(r'$-t[GeV^2/c]$')
# #     axs[2, 2].set_ylabel(r'$\delta\sigma/\sigma$')
# #     axs[2, 2].set_xlim(0, 2)
# #     axs[2, 2].set_ylim(0, 0.2)
# #     axs[2, 2].set_title('Total uncertainty')

# phi_d_2H_dsdt_results_systerr = np.zeros_like(phi_d_2H_dsdt_results_statserr[0])
# phi_d_2H_dsdt_results_totalerr = np.zeros_like(phi_d_2H_dsdt_results_statserr[0])

# for j in range(len(dsdt_index)-1):
#     phi_d_2H_dsdt_results_systerr[dsdt_index[j]:dsdt_index[j+1]] = total_uncertainties[j]
#     phi_d_2H_dsdt_results_systerr[dsdt_index[j]:dsdt_index[j+1]] = phi_d_2H_dsdt_results_systerr[dsdt_index[j]:dsdt_index[j+1]]**2

# phi_d_2H_dsdt_results_systerr += 0.0656**2  # Adding tracking efficiency
# phi_d_2H_dsdt_results_systerr += 0.05**2  # Adding photon flux uncertainty
# phi_d_2H_dsdt_results_systerr += 0.002**2  # Adding target length
# phi_d_2H_dsdt_results_systerr += 0.005**2  # Adding target density
# phi_d_2H_dsdt_results_systerr += 0.01**2  # Adding braching ratio
# phi_d_2H_dsdt_results_systerr = np.sqrt(phi_d_2H_dsdt_results_systerr)

# for j in range(len(dsdt_index)-1):
#     phi_d_2H_dsdt_results_systerr[dsdt_index[j]:dsdt_index[j+1]] = phi_d_2H_dsdt_results_systerr[dsdt_index[j]:dsdt_index[j+1]]*phi_d_2H_dsdt_results[0,dsdt_index[j]:dsdt_index[j+1]]
#     phi_d_2H_dsdt_results_totalerr[dsdt_index[j]:dsdt_index[j+1]] = np.sqrt(phi_d_2H_dsdt_results_statserr[0,dsdt_index[j]:dsdt_index[j+1]]**2 + phi_d_2H_dsdt_results_systerr[dsdt_index[j]:dsdt_index[j+1]]**2)
#     plt.scatter(phi_d_2H_dsdt_minust_center[dsdt_index[j]:dsdt_index[j+1]], phi_d_2H_dsdt_results_systerr[dsdt_index[j]:dsdt_index[j+1]]/phi_d_2H_dsdt_results[0,dsdt_index[j]:dsdt_index[j+1]], label=legend_list[j])
#     plt.legend()
#     plt.xlabel(r'$-t[GeV^2/c]$')
#     plt.ylabel(r'$\delta\sigma/\sigma$')
#     plt.xlim(0, 2)
#     plt.ylim(0, 0.2)
#     plt.title('Total point-to-point uncertainty')
# file_pdf.savefig()
# plt.close()

# fig = plt.figure(figsize=(8, 6), dpi=300)
# color_code = ['b', 'k', 'r']
# plt.errorbar(phi_d_2H_dsdt_minust_center[dsdt_index[0]:dsdt_index[1]],       phi_d_2H_dsdt_results[0, dsdt_index[0]:dsdt_index[1]],          xerr=phi_d_2H_dsdt_minust_width[dsdt_index[0]:dsdt_index[1]],        yerr=phi_d_2H_dsdt_results_statserr[0, dsdt_index[0]:dsdt_index[1]],            fmt='b.', capsize=2, capthick=1, label='This work (6-8 GeV)')
# plt.errorbar(phi_d_2H_dsdt_minust_center[dsdt_index[1]:dsdt_index[2]],       phi_d_2H_dsdt_results[0, dsdt_index[1]:dsdt_index[2]],          xerr=phi_d_2H_dsdt_minust_width[dsdt_index[1]:dsdt_index[2]],        yerr=phi_d_2H_dsdt_results_statserr[0, dsdt_index[1]:dsdt_index[2]],            fmt='k.', capsize=2, capthick=1, label='This work (8-9 GeV)')
# plt.errorbar(phi_d_2H_dsdt_minust_center[dsdt_index[2]:dsdt_index[3]],       phi_d_2H_dsdt_results[0, dsdt_index[2]:dsdt_index[3]],          xerr=phi_d_2H_dsdt_minust_width[dsdt_index[2]:dsdt_index[3]],        yerr=phi_d_2H_dsdt_results_statserr[0, dsdt_index[2]:dsdt_index[3]],            fmt='r.', capsize=2, capthick=1, label='This work (9-11 GeV)')
# plt.errorbar(phi_d_2H_dsdt_minust_center[dsdt_index[0]:dsdt_index[1]],       phi_d_2H_dsdt_results[0, dsdt_index[0]:dsdt_index[1]],          xerr=phi_d_2H_dsdt_minust_width[dsdt_index[0]:dsdt_index[1]],        yerr=phi_d_2H_dsdt_results_totalerr[dsdt_index[0]:dsdt_index[1]],            fmt='b.')
# plt.errorbar(phi_d_2H_dsdt_minust_center[dsdt_index[1]:dsdt_index[2]],       phi_d_2H_dsdt_results[0, dsdt_index[1]:dsdt_index[2]],          xerr=phi_d_2H_dsdt_minust_width[dsdt_index[1]:dsdt_index[2]],        yerr=phi_d_2H_dsdt_results_totalerr[dsdt_index[1]:dsdt_index[2]],            fmt='k.')
# plt.errorbar(phi_d_2H_dsdt_minust_center[dsdt_index[2]:dsdt_index[3]],       phi_d_2H_dsdt_results[0, dsdt_index[2]:dsdt_index[3]],          xerr=phi_d_2H_dsdt_minust_width[dsdt_index[2]:dsdt_index[3]],        yerr=phi_d_2H_dsdt_results_totalerr[dsdt_index[2]:dsdt_index[3]],            fmt='r.')
# plt.text(0.3, 0.15, 'preliminary', fontsize=15, color='r', style='italic', ha='center', va='center')
# plt.title(r"$d(\gamma, \phi d')$ differential cross section vs $-t$")
# plt.xlabel(r'$-t\ [GeV^2/c]$')
# plt.ylabel(r'$d\sigma/dt\ [nb/(GeV^2/c)]$')
# plt.xlim(0, 2)
# plt.ylim(1e-1, 1e3)
# plt.yscale('log')
# plt.legend()
# file_pdf.savefig()
# plt.close()

# # fig = plt.figure(figsize=(8, 6), dpi=300)
# # color_code = ['b', 'k', 'r']
# # # plt.errorbar(phi_d_2H_dsdt_minust_center[dsdt_index[0]:dsdt_index[1]],       phi_d_2H_dsdt_results[dsdt_index[0]:dsdt_index[1]],          xerr=phi_d_2H_dsdt_minust_width[dsdt_index[0]:dsdt_index[1]],        yerr=phi_d_2H_dsdt_results_statserr[dsdt_index[0]:dsdt_index[1]],            fmt='b.', label='This work (6-8 GeV)')
# # plt.errorbar(phi_d_2H_dsdt_minust_center[dsdt_index[1]:dsdt_index[2]],       phi_d_2H_dsdt_results[0,dsdt_index[1]:dsdt_index[2]],          xerr=phi_d_2H_dsdt_minust_width[dsdt_index[1]:dsdt_index[2]],        yerr=phi_d_2H_dsdt_results_statserr[0,dsdt_index[1]:dsdt_index[2]],            fmt='k.', label='This work (8-9 GeV)')
# # # plt.errorbar(phi_d_2H_dsdt_minust_center[dsdt_index[2]:dsdt_index[3]],       phi_d_2H_dsdt_results[dsdt_index[2]:dsdt_index[3]],          xerr=phi_d_2H_dsdt_minust_width[dsdt_index[2]:dsdt_index[3]],        yerr=phi_d_2H_dsdt_results_statserr[dsdt_index[2]:dsdt_index[3]],            fmt='r.', label='This work (9-11 GeV)')
# # plt.errorbar(np.linspace(0.1,2,191), np.loadtxt('/work/halld2/home/boyu/src_analysis/plot/vm_d/theory/output/temp/temp_sigma20_b6.5_2.22.txt')[:,2], fmt='-', color = 'y', label=r'$\sigma_{\phi N}=$20 mb, $b_{\phi N}=6.5 \rm \ GeV^{-2}, \chi^2/NDF=2.22$')
# # plt.errorbar(np.linspace(0.1,2,191), np.loadtxt('/work/halld2/home/boyu/src_analysis/plot/vm_d/theory/output/temp/temp_sigma25_b8.5_1.61.txt')[:,2], fmt='-', color = 'g', label=r'$\sigma_{\phi N}=$25 mb, $b_{\phi N}=8.5 \rm \ GeV^{-2}, \chi^2/NDF=1.61$')
# # plt.errorbar(np.linspace(0.1,2,191), np.loadtxt('/work/halld2/home/boyu/src_analysis/plot/vm_d/theory/output/temp/temp_sigma30_b10.0_1.44.txt')[:,2], fmt='-', color = 'r', label=r'$\sigma_{\phi N}=$30 mb, $b_{\phi N}=10.0 \rm \ GeV^{-2}, \chi^2/NDF=1.44$')
# # plt.errorbar(np.linspace(0.1,2,191), np.loadtxt('/work/halld2/home/boyu/src_analysis/plot/vm_d/theory/output/temp/temp_sigma35_b11.5_1.51.txt')[:,2], fmt='-', color = 'b', label=r'$\sigma_{\phi N}=$35 mb, $b_{\phi N}=11.5 \rm \ GeV^{-2}, \chi^2/NDF=1.51$')
# # plt.errorbar(np.linspace(0.1,2,191), np.loadtxt('/work/halld2/home/boyu/src_analysis/plot/vm_d/theory/output/temp/temp_sigma40_b13.0_1.69.txt')[:,2], fmt='-', color = 'c', label=r'$\sigma_{\phi N}=$40 mb, $b_{\phi N}=13.0 \rm \ GeV^{-2}, \chi^2/NDF=1.69$')
# # plt.errorbar(np.linspace(0.1,2,191), np.loadtxt('/work/halld2/home/boyu/src_analysis/plot/vm_d/theory/output/temp/temp_sigma45_b14.5_1.95.txt')[:,2], fmt='-', color = 'm', label=r'$\sigma_{\phi N}=$45 mb, $b_{\phi N}=14.5 \rm \ GeV^{-2}, \chi^2/NDF=1.95$')
# # plt.fill_between(np.linspace(0.1,2,191), np.loadtxt('/work/halld2/home/boyu/src_analysis/plot/vm_d/theory/output/temp/temp_sigma30_b10.0_1.44.txt')[:,2]-0.2*(np.loadtxt('/work/halld2/home/boyu/src_analysis/plot/vm_d/theory/output/temp/temp_sigma30_b10.0_1.44.txt')[:,2]), np.loadtxt('/work/halld2/home/boyu/src_analysis/plot/vm_d/theory/output/temp/temp_sigma30_b10.0_1.44.txt')[:,2]+0.2*(np.loadtxt('/work/halld2/home/boyu/src_analysis/plot/vm_d/theory/output/temp/temp_sigma30_b10.0_1.44.txt')[:,2]), color='r', alpha=0.2)
# # plt.text(0.3, 0.15, 'preliminary', fontsize=15, color='r', style='italic', ha='center', va='center')
# # plt.title(r"$d(\gamma, \phi d')$ differential cross section vs $-t$")
# # plt.xlabel(r'$-t\ [GeV^2/c]$')
# # plt.ylabel(r'$d\sigma/dt\ [nb/(GeV^2/c)]$')
# # plt.xlim(0, 2)
# # plt.ylim(1e-1, 1e3)
# # plt.yscale('log')
# # plt.legend()
# # file_pdf.savefig()
# # plt.close()

# # fig = plt.figure(figsize=(8, 6), dpi=300)
# # color_code = ['b', 'k', 'r']
# # plt.errorbar(phi_d_2H_dsdt_minust_center[dsdt_index[1]:dsdt_index[2]],       phi_d_2H_dsdt_results[0,dsdt_index[1]:dsdt_index[2]],          xerr=phi_d_2H_dsdt_minust_width[dsdt_index[1]:dsdt_index[2]],        yerr=phi_d_2H_dsdt_results_statserr[0,dsdt_index[1]:dsdt_index[2]],            fmt='k.', label='This work (8-9 GeV)')
# # plt.errorbar(np.linspace(0.1,2,191), np.loadtxt('/work/halld2/home/boyu/src_analysis/plot/vm_d/theory/output/case3_6.0/E_8.5_sigma_30_b_10.txt')[:,2], fmt='-', color = 'y', label=r'$\sigma_{\phi N}=$20 mb, $b_{\phi N}=6.5 \rm \ GeV^{-2}, \chi^2/NDF=2.22$')
# # plt.text(0.3, 0.15, 'preliminary', fontsize=15, color='r', style='italic', ha='center', va='center')
# # plt.title(r"$d(\gamma, \phi d')$ differential cross section vs $-t$")
# # plt.xlabel(r'$-t\ [GeV^2/c]$')
# # plt.ylabel(r'$d\sigma/dt\ [nb/(GeV^2/c)]$')
# # plt.xlim(0, 2)
# # plt.ylim(1e-1, 1e3)
# # plt.yscale('log')
# # plt.legend()
# # file_pdf.savefig()
# # plt.close()

# sphin_list = np.arange(0,40,0.5)
# bphin_list = np.arange(0,15,0.5)

# chi2_array = np.zeros((len(bphin_list), len(sphin_list)))
# for i,sphin in enumerate(sphin_list):
#     for j,bphin in enumerate(bphin_list):
#         print(f'Calculating chi2 for sigma={sphin}, b={bphin}')
#         # theory_results = np.loadtxt(f'/work/halld2/home/boyu/src_analysis/plot/vm_d/theory/output/case3_6.0/E_8.5_sigma_{sphin}_b_{bphin}.txt')
#         theory_results = np.loadtxt(f'/work/halld2/home/boyu/src_analysis/plot/vm_d/theory/output/case3_6.0/E_8.5_sigma_%.1f_b_%.1f.txt' % (sphin, bphin))
#         nuisance_opt = 0
#         for nuisance in np.arange(-10, 10.1, 0.01):
#             chi2 = 0
#             ndf = 0
#             for k in range(index[1], index[2]):
#                 t_val = phi_d_2H_dsdt_minust_center[k]
#                 data_val = phi_d_2H_dsdt_results[0,k]
#                 data_err = phi_d_2H_dsdt_results_statserr[0,k]
#                 # Find the closest theory point
#                 theory_idx = (np.abs(theory_results[:,0] - t_val)).argmin()
#                 theory_val = theory_results[theory_idx, 2]
#                 if data_err > 0 and data_val > 0 and theory_val > 0:
#                     chi2 += ((data_val - (1 + nuisance*0.2)*theory_val)**2)/(data_err**2)
#                     ndf += 1
#             chi2 += nuisance**2  # pull term
#             ndf -= 2  # two fit parameters
#             if chi2/ndf < chi2_array[j, i] or chi2_array[j, i] == 0:
#                 chi2_array[j, i] = chi2/ndf
#                 nuisance_opt = nuisance
#         # chi2_array[j, i] = np.log(chi2_array[j, i])
#         if np.abs(nuisance_opt) >= 9.9:
#             print(f'Warning: Nuisance parameter at limit for sigma={sphin}, b={bphin}, nuisance={nuisance_opt}')

# fig = plt.figure(figsize=(8, 6), dpi=300)
# # plt.contourf(sphin_list, bphin_list, chi2_array, levels=50, cmap='viridis')
# # cbar = plt.colorbar()
# # cbar.set_label(r'$\chi^2/NDF$')
# confidence_levels = [2.30, 6.18, 11.83]  # 68%, 95%, 99.7% for 2 parameters
# best_fit_idx = np.unravel_index(np.argmin(chi2_array), chi2_array.shape)
# best_chi2 = chi2_array[best_fit_idx]
# for level in confidence_levels:
#     contour_level = best_chi2 + level/ndf
#     plt.contour(sphin_list, bphin_list, chi2_array, levels=[contour_level], colors='black', linestyles='dashed')
#     # plt.text(best_fit_idx[1]*5 + 5, best_fit_idx[0]*5 + 20 + level, f'{int(level*100)/100}', color='black')
# plt.xlabel(r'$\sigma_{\phi N}\ [mb]$')
# plt.ylabel(r'$\rm b_{\phi N}\ [GeV^{-2}]$')
# plt.title(r'$\chi^2/NDF$ map for $\phi-D$ scattering parameters')
# file_pdf.savefig()
# plt.close()

# # sphin_list = np.arange(0,120,1)
# # bphin_list = np.arange(0,30,1)

# # chi2_array = np.zeros((len(bphin_list), len(sphin_list)))
# # for i,sphin in enumerate(sphin_list):
# #     for j,bphin in enumerate(bphin_list):
# #         theory_results = np.loadtxt(f'/work/halld2/home/boyu/src_analysis/plot/vm_d/theory/output/case3_4.0_alpha0_61237/E_2.1_sigma_%.0f_b_%.0f.txt' % (sphin, bphin))
# #         chi2 = 0
# #         ndf = 0
# #         for k in range(len(phi_d_2H_dsdt_clas_minust_center)):
# #             t_val = phi_d_2H_dsdt_clas_minust_center[k]
# #             data_val = phi_d_2H_dsdt_clas_results_16[k]
# #             data_err = phi_d_2H_dsdt_clas_results_16_totalerr[k]
# #             # Find the closest theory point
# #             theory_idx = (np.abs(theory_results[:,0] - t_val)).argmin()
# #             theory_val = theory_results[theory_idx, 2]
# #             if data_err > 0 and data_val > 0 and theory_val > 0:
# #                 chi2 += ((data_val - theory_val)**2)/(data_err**2 + (0.1*theory_val)**2)  # adding 20% theory uncertainty
# #                 ndf += 1
# #         # ndf -= 2  # two fit parameters
# #         chi2_array[j, i] = np.log(chi2/ndf)

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


# # fig = plt.figure(figsize=(8, 6), dpi=300)
# # color_code = ['b', 'k', 'r']
# # plt.errorbar(phi_d_2H_dsdt_minust_center[index[0]:index[1]],       phi_d_2H_dsdt_results[index[0]:index[1]],          xerr=phi_d_2H_dsdt_minust_width[index[0]:index[1]],        yerr=phi_d_2H_dsdt_results_statserr[index[0]:index[1]],            fmt='b.', label='This work (6-8 GeV)')
# # plt.errorbar(phi_d_2H_dsdt_minust_center[index[1]:index[2]],       phi_d_2H_dsdt_results[index[1]:index[2]],          xerr=phi_d_2H_dsdt_minust_width[index[1]:index[2]],        yerr=phi_d_2H_dsdt_results_statserr[index[1]:index[2]],            fmt='k.', label='This work (8-9 GeV)')
# # plt.errorbar(phi_d_2H_dsdt_minust_center[index[2]:index[3]],       phi_d_2H_dsdt_results[index[2]:index[3]],          xerr=phi_d_2H_dsdt_minust_width[index[2]:index[3]],        yerr=phi_d_2H_dsdt_results_statserr[index[2]:index[3]],            fmt='r.', label='This work (9-11 GeV)')
# # plt.errorbar(np.linspace(0.1,2,191), np.loadtxt('/work/halld2/home/boyu/src_analysis/plot/vm_d/theory/8.5gev_sigma10_b6_case3_new.txt')[:,2], fmt='-', color = 'g', label='E=8.5 GeV $\sigma_{\phi N}=$10 mb (VMD)')
# # plt.fill_between(np.linspace(0.1,2,191), np.loadtxt('/work/halld2/home/boyu/src_analysis/plot/vm_d/theory/8.5gev_sigma10_b6_case3_new.txt')[:,2]-0.2*(np.loadtxt('/work/halld2/home/boyu/src_analysis/plot/vm_d/theory/8.5gev_sigma10_b6_case3_new.txt')[:,2]), np.loadtxt('/work/halld2/home/boyu/src_analysis/plot/vm_d/theory/8.5gev_sigma10_b6_case3_new.txt')[:,2]+0.2*(np.loadtxt('/work/halld2/home/boyu/src_analysis/plot/vm_d/theory/8.5gev_sigma10_b6_case3_new.txt')[:,2]), color='g', alpha=0.2)
# # plt.errorbar(np.linspace(0.1,2,191), np.loadtxt('/work/halld2/home/boyu/src_analysis/plot/vm_d/theory/8.5gev_sigma30_b8_case3.txt')[:,2], fmt='-', color = 'y', label=r'E=8.5 GeV $\sigma_{\phi N}=$30 mb ($b_{\phi N}=8 \rm \ GeV^{-2}$)')
# # plt.fill_between(np.linspace(0.1,2,191), np.loadtxt('/work/halld2/home/boyu/src_analysis/plot/vm_d/theory/8.5gev_sigma30_b8_case3.txt')[:,2]-0.2*(np.loadtxt('/work/halld2/home/boyu/src_analysis/plot/vm_d/theory/8.5gev_sigma30_b8_case3.txt')[:,2]), np.loadtxt('/work/halld2/home/boyu/src_analysis/plot/vm_d/theory/8.5gev_sigma30_b8_case3.txt')[:,2]+0.2*(np.loadtxt('/work/halld2/home/boyu/src_analysis/plot/vm_d/theory/8.5gev_sigma30_b8_case3.txt')[:,2]), color='y', alpha=0.2)
# # plt.text(0.3, 0.15, 'preliminary', fontsize=15, color='r', style='italic', ha='center', va='center')
# # plt.title(r"$d(\gamma, \phi d')$ differential cross section vs $-t$")
# # plt.xlabel(r'$-t\ [GeV^2/c]$')
# # plt.ylabel(r'$d\sigma/dt\ [nb/(GeV^2/c)]$')
# # plt.xlim(0, 2)
# # plt.ylim(1e-1, 1e3)
# # plt.yscale('log')
# # plt.legend()
# # file_pdf.savefig()
# # plt.close()

# # fig = plt.figure(figsize=(8, 6), dpi=300)
# # color_code = ['b', 'k', 'r']
# # plt.errorbar(phi_d_2H_dsdt_minust_center[index[0]:index[1]],       phi_d_2H_dsdt_results[index[0]:index[1]],          xerr=phi_d_2H_dsdt_minust_width[index[0]:index[1]],        yerr=phi_d_2H_dsdt_results_statserr[index[0]:index[1]],            fmt='b.', label='6-8 GeV')
# # plt.errorbar(phi_d_2H_dsdt_minust_center[index[1]:index[2]],       phi_d_2H_dsdt_results[index[1]:index[2]],          xerr=phi_d_2H_dsdt_minust_width[index[1]:index[2]],        yerr=phi_d_2H_dsdt_results_statserr[index[1]:index[2]],            fmt='k.', label='8-9 GeV')
# # plt.errorbar(phi_d_2H_dsdt_minust_center[index[2]:index[3]],       phi_d_2H_dsdt_results[index[2]:index[3]],          xerr=phi_d_2H_dsdt_minust_width[index[2]:index[3]],        yerr=phi_d_2H_dsdt_results_statserr[index[2]:index[3]],            fmt='r.', label='9-11 GeV')
# # plt.plot(np.linspace(0.1,2,191), np.loadtxt('/work/halld2/home/boyu/src_analysis/plot/vm_d/theory/8.5gev_sigma30_b4_case3_new.txt')[:,2], '-', color = 'k', label='Theory 8.5 GeV 30mb b=4')
# # plt.errorbar(np.linspace(0.1,2,191), np.loadtxt('/work/halld2/home/boyu/src_analysis/plot/vm_d/theory/8.5gev_sigma30_b6_case3_new.txt')[:,2], yerr=0.2*(np.loadtxt('/work/halld2/home/boyu/src_analysis/plot/vm_d/theory/8.5gev_sigma30_b6_case3_new.txt')[:,2]), fmt='-', color = 'b', label='Theory 8.5 GeV 30mb b=6')
# # plt.plot(np.linspace(0.1,2,191), np.loadtxt('/work/halld2/home/boyu/src_analysis/plot/vm_d/theory/8.5gev_sigma30_b8_case3_new.txt')[:,2], '-', color = 'r', label='Theory 8.5 GeV 30mb b=8')
# # plt.plot(np.linspace(0.1,2,191), np.loadtxt('/work/halld2/home/boyu/src_analysis/plot/vm_d/theory/8.5gev_sigma30_b4_10_case3.txt')[:,2], '-', color = 'm', label='Theory 8.5 GeV 30mb b=4_10')
# # plt.plot(np.linspace(0.1,2,191), np.loadtxt('/work/halld2/home/boyu/src_analysis/plot/vm_d/theory/8.5gev_sigma30_b10_case3_new.txt')[:,2], '-', color = 'g', label='Theory 8.5 GeV 30mb b=10')
# # plt.text(0.3, 0.15, 'preliminary', fontsize=15, color='r', style='italic', ha='center', va='center')
# # plt.title(r"$d(\gamma, \phi d')$ differential cross section vs $-t$")
# # plt.xlabel(r'$-t\ [GeV^2/c]$')
# # plt.ylabel(r'$d\sigma/dt\ [nb/(GeV^2/c)]$')
# # plt.xlim(0, 2)
# # plt.ylim(1e-1, 1e3)
# # plt.yscale('log')
# # plt.legend()
# # file_pdf.savefig()
# # plt.close()

# fig = plt.figure(figsize=(8, 6), dpi=300)
# color_code = ['b', 'k', 'r']
# plt.errorbar(phi_d_2H_dsdt_minust_center[index[0]:index[1]],       phi_d_2H_dsdt_results[0, index[0]:index[1]],          xerr=phi_d_2H_dsdt_minust_width[index[0]:index[1]],        yerr=phi_d_2H_dsdt_results_statserr[0, index[0]:index[1]],            fmt='b.', label='6-8 GeV')
# plt.errorbar(phi_d_2H_dsdt_minust_center[index[1]:index[2]],       phi_d_2H_dsdt_results[0, index[1]:index[2]],          xerr=phi_d_2H_dsdt_minust_width[index[1]:index[2]],        yerr=phi_d_2H_dsdt_results_statserr[0, index[1]:index[2]],            fmt='k.', label='8-9 GeV')
# plt.errorbar(phi_d_2H_dsdt_minust_center[index[2]:index[3]],       phi_d_2H_dsdt_results[0, index[2]:index[3]],          xerr=phi_d_2H_dsdt_minust_width[index[2]:index[3]],        yerr=phi_d_2H_dsdt_results_statserr[0, index[2]:index[3]],            fmt='r.', label='9-11 GeV')
# for i in range(3):
#     curve_fit_params, curve_fit_cov = curve_fit(dsdt_func, phi_d_2H_dsdt_minust_center[index[i]:index[i+1]], phi_d_2H_dsdt_results[0, index[i]:index[i+1]], p0=[3000, 15, 20, 3])
#     curve_fit_residuals = phi_d_2H_dsdt_results[0, index[i]:index[i+1]] - dsdt_func(phi_d_2H_dsdt_minust_center[index[i]:index[i+1]], curve_fit_params[0], curve_fit_params[1], curve_fit_params[2], curve_fit_params[3])
#     reduced_chi2 = np.sum((curve_fit_residuals/phi_d_2H_dsdt_results_statserr[0, index[i]:index[i+1]])**2)/(len(phi_d_2H_dsdt_results[0, index[i]:index[i+1]])-4)
#     plt.plot(np.linspace(0, 2, 100), dsdt_func(np.linspace(0, 2, 100), curve_fit_params[0], curve_fit_params[1], curve_fit_params[2], curve_fit_params[3]), '--', color = color_code[i], label='Paras: %.2f, %.2f, %.2f, %.2f, $\chi^2$/ndf = %.2f' % (curve_fit_params[0], curve_fit_params[1], curve_fit_params[2], curve_fit_params[3], reduced_chi2))
# plt.text(0.3, 0.15, 'preliminary', fontsize=15, color='r', style='italic', ha='center', va='center')
# plt.title(r"$d(\gamma, \phi d')$ differential cross section vs $-t$")
# plt.xlabel(r'$-t\ [GeV^2/c]$')
# plt.ylabel(r'$d\sigma/dt\ [nb/(GeV^2/c)]$')
# plt.xlim(0, 2)
# plt.ylim(1e-1, 1e3)
# plt.yscale('log')
# plt.legend()
# file_pdf.savefig()
# plt.close()

# fig = plt.figure(figsize=(8, 6), dpi=300)
# color_code = ['b', 'k', 'r']
# plt.errorbar(phi_d_2H_dsdt_minust_center[index[0]:index[1]],       phi_d_2H_dsdt_results[index[0]:index[1]],          xerr=phi_d_2H_dsdt_minust_width[index[0]:index[1]],        yerr=phi_d_2H_dsdt_results_statserr[index[0]:index[1]],            fmt='b.', label='6-8 GeV')
# plt.errorbar(phi_d_2H_dsdt_minust_center[index[1]:index[2]],       10*phi_d_2H_dsdt_results[index[1]:index[2]],          xerr=phi_d_2H_dsdt_minust_width[index[1]:index[2]],        yerr=10*phi_d_2H_dsdt_results_statserr[index[1]:index[2]],            fmt='k.', label='8-9 GeV')
# plt.errorbar(phi_d_2H_dsdt_minust_center[index[2]:index[3]],       100*phi_d_2H_dsdt_results[index[2]:index[3]],          xerr=phi_d_2H_dsdt_minust_width[index[2]:index[3]],        yerr=100*phi_d_2H_dsdt_results_statserr[index[2]:index[3]],            fmt='r.', label='9-11 GeV')
# for i in range(3):
#     curve_fit_params, curve_fit_cov = curve_fit(dsdt_func, phi_d_2H_dsdt_minust_center[index[i]:index[i+1]], phi_d_2H_dsdt_results[index[i]:index[i+1]], p0=[0, 0, 0, 0])
#     curve_fit_residuals = phi_d_2H_dsdt_results[index[i]:index[i+1]] - dsdt_func(phi_d_2H_dsdt_minust_center[index[i]:index[i+1]], curve_fit_params[0], curve_fit_params[1], curve_fit_params[2], curve_fit_params[3])
#     reduced_chi2 = np.sum((curve_fit_residuals/phi_d_2H_dsdt_results_statserr[index[i]:index[i+1]])**2)/(len(phi_d_2H_dsdt_results[index[i]:index[i+1]])-4)
#     plt.plot(np.linspace(0, 2, 100), pow(10,i)*dsdt_func(np.linspace(0, 2, 100), curve_fit_params[0], curve_fit_params[1], curve_fit_params[2], curve_fit_params[3]), '--', color = color_code[i])
# plt.text(0.3, 0.15, 'preliminary', fontsize=15, color='r', style='italic', ha='center', va='center')
# plt.title(r"$d(\gamma, \phi d')$ differential cross section vs $-t$")
# plt.xlabel(r'$-t\ [GeV^2/c]$')
# plt.ylabel(r'$d\sigma/dt\ [nb/(GeV^2/c)]$')
# plt.xlim(0, 2)
# plt.ylim(1e-1, 1e4)
# plt.yscale('log')
# plt.legend()
# file_pdf.savefig()
# plt.close()

# fig = plt.figure(figsize=(8, 6), dpi=300)
# plt.errorbar(phi_d_2H_dsdt_clas_minust_center,  phi_d_2H_dsdt_clas_results_16,  yerr=phi_d_2H_dsdt_clas_results_16_statserr,    fmt='yo', markersize=3, label='CLAS 1.6-2.6 GeV')
# plt.errorbar(phi_d_2H_dsdt_clas_minust_center,  phi_d_2H_dsdt_clas_results_26,  yerr=phi_d_2H_dsdt_clas_results_26_statserr,    fmt='co', markersize=3, label='CLAS 2.6-3.6 GeV')
# curve_fit_params, curve_fit_cov = curve_fit(dsdt_func, phi_d_2H_dsdt_clas_minust_center, phi_d_2H_dsdt_clas_results_16, p0=[0, 0, 0, 0])
# curve_fit_residuals = phi_d_2H_dsdt_clas_results_16 - dsdt_func(phi_d_2H_dsdt_clas_minust_center, curve_fit_params[0], curve_fit_params[1], curve_fit_params[2], curve_fit_params[3])
# reduced_chi2 = np.sum((curve_fit_residuals/phi_d_2H_dsdt_clas_results_16_statserr)**2)/(len(phi_d_2H_dsdt_clas_results_16)-4)
# print(curve_fit_params)
# plt.plot(np.linspace(0, 2, 100), dsdt_func(np.linspace(0, 2, 100), curve_fit_params[0], curve_fit_params[1], curve_fit_params[2], curve_fit_params[3]), '--', color = 'y')
# curve_fit_params, curve_fit_cov = curve_fit(dsdt_func, phi_d_2H_dsdt_clas_minust_center, phi_d_2H_dsdt_clas_results_26, p0=[0, 2, 0, 8], bounds=([-np.inf, 0, -np.inf, 0], [np.inf, 5, np.inf, 10]))
# curve_fit_residuals = phi_d_2H_dsdt_clas_results_26 - dsdt_func(phi_d_2H_dsdt_clas_minust_center, curve_fit_params[0], curve_fit_params[1], curve_fit_params[2], curve_fit_params[3])
# reduced_chi2 = np.sum((curve_fit_residuals/phi_d_2H_dsdt_clas_results_26_statserr)**2)/(len(phi_d_2H_dsdt_clas_results_26)-4)
# plt.plot(np.linspace(0, 2, 100), dsdt_func(np.linspace(0, 2, 100), curve_fit_params[0], curve_fit_params[1], curve_fit_params[2], curve_fit_params[3]), '--', color = 'c')
# print(curve_fit_params)
# plt.title(r"$d(\gamma, \phi d')$ differential cross section vs $-t$")
# plt.xlabel(r'$-t\ [GeV^2/c]$')
# plt.ylabel(r'$d\sigma/dt\ [nb/(GeV^2/c)]$')
# plt.xlim(0, 2)
# plt.ylim(1e-1, 1e4)
# plt.yscale('log')
# plt.legend()
# file_pdf.savefig()
# plt.close()

# iteration_index = 6

# # Read the bin edges
# phi_d_2H_dsdt_energy_low_simweight          = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_dsdt_nominal.txt')[:,2]
# phi_d_2H_dsdt_energy_high_simweight         = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_dsdt_nominal.txt')[:,3]
# phi_d_2H_dsdt_energy_center_simweight       = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_dsdt_nominal.txt')[:,0]
# phi_d_2H_dsdt_energy_width_simweight        = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_dsdt_nominal.txt')[:,1]
# phi_d_2H_dsdt_minust_low_simweight          = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_dsdt_nominal.txt')[:,6]
# phi_d_2H_dsdt_minust_high_simweight         = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_dsdt_nominal.txt')[:,7]
# phi_d_2H_dsdt_minust_center_simweight       = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_dsdt_nominal.txt')[:,4]
# phi_d_2H_dsdt_minust_width_simweight        = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_dsdt_nominal.txt')[:,5]

# # Read the yield numbers
# phi_d_2H_dsdt_yield_data_simweight            = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_dsdt_nominal.txt')[:,8]
# phi_d_2H_dsdt_yield_data_statserr_simweight   = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_dsdt_nominal.txt')[:,9]/phi_d_2H_dsdt_yield_data_simweight
# phi_d_2H_dsdt_yield_sim_simweight             = np.loadtxt(f'output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_dsdt_simweight_pass{iteration_index}.txt')[:,8]
# phi_d_2H_dsdt_yield_sim_statser_simweight     = np.loadtxt(f'output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_dsdt_simweight_pass{iteration_index}.txt')[:,9]/phi_d_2H_dsdt_yield_sim_simweight
# phi_d_2H_dsdt_yield_tagged_simweight          = np.loadtxt(f'output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_dsdt_simweight_pass{iteration_index}.txt')[:,8]
# phi_d_2H_dsdt_yield_tagged_statserr_simweight = np.loadtxt(f'output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_dsdt_simweight_pass{iteration_index}.txt')[:,9]/phi_d_2H_dsdt_yield_tagged_simweight
# phi_d_2H_dsdt_efficiency_simweight            = phi_d_2H_dsdt_yield_sim_simweight/phi_d_2H_dsdt_yield_tagged_simweight
# phi_d_2H_dsdt_efficiency_statserr_simweight   = np.sqrt(phi_d_2H_dsdt_yield_sim_statser_simweight**2 + phi_d_2H_dsdt_yield_tagged_statserr_simweight**2)
# phi_d_2H_dsdt_results_simweight               = phi_d_2H_dsdt_yield_data_simweight/phi_d_2H_dsdt_efficiency_simweight/lumi(phi_d_2H_dsdt_energy_low_simweight, phi_d_2H_dsdt_energy_high_simweight, 28)/(phi_d_2H_dsdt_minust_high_simweight-phi_d_2H_dsdt_minust_low_simweight)/0.489/1000
# phi_d_2H_dsdt_results_statserr_simweight      = phi_d_2H_dsdt_results_simweight*np.sqrt(phi_d_2H_dsdt_yield_data_statserr_simweight**2 + phi_d_2H_dsdt_efficiency_statserr_simweight**2)

# # Find the indices for the different energy and t bins
# index = []
# for i in range(len(phi_d_2H_dsdt_results_simweight)):
#     if (i == 0):
#         index.append(i)
#     elif (i == len(phi_d_2H_dsdt_results_simweight) - 1):
#         index.append(i+1)
#     else:
#         if (phi_d_2H_dsdt_energy_low_simweight[i] != phi_d_2H_dsdt_energy_low_simweight[i-1]):
#             index.append(i)

# fig = plt.figure(figsize=(8, 6), dpi=300)
# color_code = ['b', 'k', 'r']
# plt.errorbar(phi_d_2H_dsdt_minust_center_simweight[index[0]:index[1]],       phi_d_2H_dsdt_results_simweight[index[0]:index[1]],          xerr=phi_d_2H_dsdt_minust_width_simweight[index[0]:index[1]],        yerr=phi_d_2H_dsdt_results_statserr_simweight[index[0]:index[1]],            fmt='b.', label='6-8 GeV')
# plt.errorbar(phi_d_2H_dsdt_minust_center_simweight[index[1]:index[2]],       phi_d_2H_dsdt_results_simweight[index[1]:index[2]],          xerr=phi_d_2H_dsdt_minust_width_simweight[index[1]:index[2]],        yerr=phi_d_2H_dsdt_results_statserr_simweight[index[1]:index[2]],            fmt='k.', label='8-9 GeV')
# plt.errorbar(phi_d_2H_dsdt_minust_center_simweight[index[2]:index[3]],       phi_d_2H_dsdt_results_simweight[index[2]:index[3]],          xerr=phi_d_2H_dsdt_minust_width_simweight[index[2]:index[3]],        yerr=phi_d_2H_dsdt_results_statserr_simweight[index[2]:index[3]],            fmt='r.', label='9-11 GeV')
# for i in range(3):
#     curve_fit_params, curve_fit_cov = curve_fit(dsdt_func, phi_d_2H_dsdt_minust_center_simweight[index[i]:index[i+1]], phi_d_2H_dsdt_results_simweight[index[i]:index[i+1]], p0=[3000, 15, 20, 3])
#     curve_fit_residuals = phi_d_2H_dsdt_results_simweight[index[i]:index[i+1]] - dsdt_func(phi_d_2H_dsdt_minust_center_simweight[index[i]:index[i+1]], curve_fit_params[0], curve_fit_params[1], curve_fit_params[2], curve_fit_params[3])
#     reduced_chi2 = np.sum((curve_fit_residuals/phi_d_2H_dsdt_results_statserr_simweight[index[i]:index[i+1]])**2)/(len(phi_d_2H_dsdt_results_simweight[index[i]:index[i+1]])-4)
#     plt.plot(np.linspace(0, 2, 100), dsdt_func(np.linspace(0, 2, 100), curve_fit_params[0], curve_fit_params[1], curve_fit_params[2], curve_fit_params[3]), '--', color = color_code[i], label='Paras: %.3e, %.3e, %.3e, %.3e, $\chi^2$/ndf = %.2f' % (curve_fit_params[0], curve_fit_params[1], curve_fit_params[2], curve_fit_params[3], reduced_chi2))
#     print('Energy bin %d: a1 = %.3e +/- %.3e, b1 = %.3e +/- %.3e, a2 = %.3e +/- %.3e, b2 = %.3e +/- %.3e, chi2/ndf = %.2f' % (i, curve_fit_params[0], np.sqrt(curve_fit_cov[0][0]), curve_fit_params[1], np.sqrt(curve_fit_cov[1][1]), curve_fit_params[2], np.sqrt(curve_fit_cov[2][2]), curve_fit_params[3], np.sqrt(curve_fit_cov[3][3]), reduced_chi2))
# plt.text(0.3, 0.15, 'preliminary', fontsize=15, color='r', style='italic', ha='center', va='center')
# plt.title(r"$d(\gamma, \phi d')$ differential cross section vs $-t$")
# plt.xlabel(r'$-t\ [GeV^2/c]$')
# plt.ylabel(r'$d\sigma/dt\ [nb/(GeV^2/c)]$')
# plt.xlim(0, 2)
# plt.ylim(1e-1, 1e3)
# plt.yscale('log')
# plt.legend()
# file_pdf.savefig()
# plt.close()


# a1_list = np.array([3632.01, 6127.78, 6417.46, 6442.44, 6444.93, 6444.97, 6445.03])
# b1_list = np.array([15.82, 17.10, 17.21, 17.22, 17.22, 17.22, 17.22])
# a2_list = np.array([12.19, 12.43, 12.42, 12.42, 12.42, 12.42, 12.42])
# b2_list = np.array([2.49, 2.50, 2.49, 2.49, 2.49, 2.49, 2.49])

# # plt.plot(np.arange(len(a1_list)), (a1_list-a1_list[-1])/a1_list[-1], 'bo-', label='a1')
# # plt.plot(np.arange(len(b1_list)), (b1_list-b1_list[-1])/b1_list[-1], 'ro-', label='b1')
# # plt.plot(np.arange(len(a2_list)), (a2_list-a2_list[-1])/a2_list[-1], 'go-', label='a2')
# # plt.plot(np.arange(len(b2_list)), (b2_list-b2_list[-1])/b2_list[-1], 'mo-', label='b2')
# plt.plot(np.arange(1, len(a1_list)), (a1_list[0:-1]-a1_list[1:])/a1_list[0:-1], 'bo-', label='a1')
# plt.plot(np.arange(1, len(b1_list)), (b1_list[0:-1]-b1_list[1:])/b1_list[0:-1], 'ro-', label='b1')
# plt.plot(np.arange(1, len(a2_list)), (a2_list[0:-1]-a2_list[1:])/a2_list[0:-1], 'go-', label='a2')
# plt.plot(np.arange(1, len(b2_list)), (b2_list[0:-1]-b2_list[1:])/b2_list[0:-1], 'mo-', label='b2')
# plt.axhline(y=0, color='k', linestyle='--')
# plt.title(r"Fit parameter variation with different simulation weights")
# plt.xlabel('Iteration')
# plt.ylabel('Relative variation')
# plt.legend()
# file_pdf.savefig()
# plt.close()

# # Data points from CLAS
# phi_d_2H_dsdt_clas_minust_low           = np.array([0.350, 0.375, 0.400, 0.425, 0.450, 0.500, 0.550, 0.600, 0.700, 0.800, 1.000, 1.200, 1.400])
# phi_d_2H_dsdt_clas_minust_high          = np.array([0.375, 0.400, 0.425, 0.450, 0.500, 0.550, 0.600, 0.700, 0.800, 1.000, 1.200, 1.400, 2.000])
# phi_d_2H_dsdt_clas_minust_middle        = (phi_d_2H_dsdt_clas_minust_high + phi_d_2H_dsdt_clas_minust_low) / 2
# phi_d_2H_dsdt_clas_minust_size          = (phi_d_2H_dsdt_clas_minust_high - phi_d_2H_dsdt_clas_minust_low) / 2
# phi_d_2H_dsdt_clas_minust_center        = np.array([0.360, 0.385, 0.410, 0.435, 0.474, 0.524, 0.574, 0.646, 0.746, 0.888, 1.091, 1.292, 1.637])
# phi_d_2H_dsdt_clas_results_16           = np.array([10.21, 8.85, 7.32, 6.16, 4.73, 3.52, 2.66, 2.17, 1.40, 0.94, 0.57, 0.28, 0.19])
# phi_d_2H_dsdt_clas_results_16_statserr  = np.array([0.82, 0.75, 0.59, 0.55, 0.34, 0.28, 0.24, 0.15, 0.12, 0.07, 0.06, 0.05, 0.02])
# phi_d_2H_dsdt_clas_results_16_systerr   = np.array([1.70, 1.11, 0.94, 0.81, 0.60, 0.51, 0.38, 0.26, 0.16, 0.11, 0.07, 0.04, 0.03])
# phi_d_2H_dsdt_clas_results_16_totalerr  = np.sqrt(phi_d_2H_dsdt_clas_results_16_statserr**2 + phi_d_2H_dsdt_clas_results_16_systerr**2)
# phi_d_2H_dsdt_clas_results_26           = np.array([8.63, 6.80, 4.57, 5.76, 3.99, 3.59, 2.11, 1.83, 1.32, 0.96, 0.57, 0.36, 0.15])
# phi_d_2H_dsdt_clas_results_26_statserr  = np.array([0.80, 0.69, 0.53, 0.56, 0.33, 0.29, 0.22, 0.14, 0.12, 0.07, 0.05, 0.04, 0.02])
# phi_d_2H_dsdt_clas_results_26_systerr   = np.array([1.04, 1.07, 0.74, 0.65, 0.55, 0.55, 0.28, 0.24, 0.20, 0.11, 0.06, 0.05, 0.02])
# phi_d_2H_dsdt_clas_results_26_totalerr  = np.sqrt(phi_d_2H_dsdt_clas_results_26_statserr**2 + phi_d_2H_dsdt_clas_results_26_systerr**2)

# # Data points from LEPS
# phi_d_2H_dsdt_leps_minust_low           = np.array([0.38, 0.36, 0.34, 0.32, 0.30, 0.28, 0.26, 0.24, 0.22, 0.20, 0.18, 0.16, 0.14, 0.12, 0.10])
# phi_d_2H_dsdt_leps_minust_high          = np.array([0.4, 0.38, 0.36, 0.34, 0.32, 0.30, 0.28, 0.26, 0.24, 0.22, 0.20, 0.18, 0.16, 0.14, 0.12])
# phi_d_2H_dsdt_leps_minust_center        = (phi_d_2H_dsdt_leps_minust_high + phi_d_2H_dsdt_leps_minust_low) / 2
# phi_d_2H_dsdt_leps_minust_width         = (phi_d_2H_dsdt_leps_minust_high - phi_d_2H_dsdt_leps_minust_low) / 2
# phi_d_2H_dsdt_leps_results_157          = np.array([0.0005, 0.004, 0.0087, 0.0068, 0.0238, 0.0317, 0.0567, 0.0722, 0.092, 0.1186, 0.1749, 0.2033, 0.2544, 0.3101, 0.3396])*1000
# phi_d_2H_dsdt_leps_results_157_statserr = np.array([0.0005, 0.002, 0.0035, 0.0028, 0.007, 0.0076, 0.0102, 0.0118, 0.0142, 0.0137, 0.0159, 0.0148, 0.0166, 0.0152, 0.0143])*1000


################################################################# OBSERVABLES #############################################################################################################################
# fig = plt.figure(figsize=(8, 6), dpi=300)
# plt.errorbar(phi_d_2H_dsdt_clas_minust_center,  phi_d_2H_dsdt_clas_results_16,  yerr=phi_d_2H_dsdt_clas_results_16_statserr,    fmt='ys', markerfacecolor='white', markersize=3, capsize=2, capthick=1)
# plt.errorbar(phi_d_2H_dsdt_clas_minust_center,  phi_d_2H_dsdt_clas_results_16,  yerr=phi_d_2H_dsdt_clas_results_16_totalerr,    fmt='ys', markerfacecolor='white', markersize=3, label='CLAS (1.6-2.6 GeV)')
# plt.errorbar(phi_d_2H_dsdt_clas_minust_center,  phi_d_2H_dsdt_clas_results_26,  yerr=phi_d_2H_dsdt_clas_results_26_statserr,    fmt='cs', markerfacecolor='white', markersize=3, capsize=2, capthick=1)
# plt.errorbar(phi_d_2H_dsdt_clas_minust_center,  phi_d_2H_dsdt_clas_results_26,  yerr=phi_d_2H_dsdt_clas_results_26_totalerr,    fmt='cs', markerfacecolor='white', markersize=3, label='CLAS (2.6-3.6 GeV)')
# plt.errorbar(phi_d_2H_dsdt_leps_minust_center,  phi_d_2H_dsdt_leps_results_157, yerr=phi_d_2H_dsdt_leps_results_157_statserr,   fmt='g^', markerfacecolor='white', markersize=3, capsize=2, capthick=1)
# plt.errorbar(phi_d_2H_dsdt_leps_minust_center,  phi_d_2H_dsdt_leps_results_157, yerr=phi_d_2H_dsdt_leps_results_157_statserr,   fmt='g^', markerfacecolor='white', markersize=3, label='LEPS (1.57-2.37 GeV)')
# plt.errorbar(phi_d_2H_dsdt_minust_center[dsdt_index[0]:dsdt_index[1]], phi_d_2H_dsdt_results[0,dsdt_index[0]:dsdt_index[1]], yerr=phi_d_2H_dsdt_results_statserr[0,dsdt_index[0]:dsdt_index[1]], fmt='bo', markersize=3, capsize=2, capthick=1)
# plt.errorbar(phi_d_2H_dsdt_minust_center[dsdt_index[0]:dsdt_index[1]], phi_d_2H_dsdt_results[0,dsdt_index[0]:dsdt_index[1]], yerr=phi_d_2H_dsdt_results_statserr[0,dsdt_index[0]:dsdt_index[1]], fmt='bo', markersize=3, label='This work (6-8 GeV)')
# plt.errorbar(phi_d_2H_dsdt_minust_center[dsdt_index[1]:dsdt_index[2]], phi_d_2H_dsdt_results[0,dsdt_index[1]:dsdt_index[2]], yerr=phi_d_2H_dsdt_results_statserr[0,dsdt_index[1]:dsdt_index[2]], fmt='ko', markersize=3, capsize=2, capthick=1)
# plt.errorbar(phi_d_2H_dsdt_minust_center[dsdt_index[1]:dsdt_index[2]], phi_d_2H_dsdt_results[0,dsdt_index[1]:dsdt_index[2]], yerr=phi_d_2H_dsdt_results_statserr[0,dsdt_index[1]:dsdt_index[2]], fmt='ko', markersize=3, label='This work (8-9 GeV)')
# plt.errorbar(phi_d_2H_dsdt_minust_center[dsdt_index[2]:dsdt_index[3]], phi_d_2H_dsdt_results[0,dsdt_index[2]:dsdt_index[3]], yerr=phi_d_2H_dsdt_results_statserr[0,dsdt_index[2]:dsdt_index[3]], fmt='ro', markersize=3, capsize=2, capthick=1)
# plt.errorbar(phi_d_2H_dsdt_minust_center[dsdt_index[2]:dsdt_index[3]], phi_d_2H_dsdt_results[0,dsdt_index[2]:dsdt_index[3]], yerr=phi_d_2H_dsdt_results_statserr[0,dsdt_index[2]:dsdt_index[3]], fmt='ro', markersize=3, label='This work (9-11 GeV)')
# plt.text(0.3, 0.15, 'preliminary', fontsize=15, color='r', style='italic', ha='center', va='center')
# plt.title(r"$d(\gamma, \phi d')$ differential cross section vs $-t$")
# plt.xlabel(r'$-t\ [GeV^2/c]$')
# plt.ylabel(r'$d\sigma/dt\ [nb/(GeV^2/c)]$')
# plt.xlim(0, 2)
# plt.ylim(1e-1, 1e3)
# plt.yscale('log')
# plt.legend()
# file_pdf.savefig()
# plt.close()


# #======================================================================phi_d_2H_Wcostheta======================================================================

# # Plot the results
# fig = plt.figure(figsize=(6*(len(index) - 1), 6), dpi=300)
# gs = fig.add_gridspec(1, len(index) - 1, wspace=0)
# axs = gs.subplots(sharex=True, sharey=True)
# for i in range(len(index) - 1):
#     axs[i].errorbar(phi_d_2H_Wcostheta_costheta_center[index[i]:index[i+1]], phi_d_2H_Wcostheta_results[index[i]:index[i+1]], xerr=phi_d_2H_Wcostheta_costheta_width[index[i]:index[i+1]], yerr=phi_d_2H_Wcostheta_results_statserr[index[i]:index[i+1]], fmt='k.', label='This work')
#     axs[i].set_title(r'$%.1f<\mathrm{E}_{\gamma}<%.1f\ \mathrm{GeV},\  %.1f<-t<%.1f\ \mathrm{GeV}^2$' % (phi_d_2H_Wcostheta_energy_low[index[i]], phi_d_2H_Wcostheta_energy_high[index[i]], phi_d_2H_Wcostheta_minust_low[index[i]], phi_d_2H_Wcostheta_minust_high[index[i]]))
#     curve_fit_params, curve_fit_cov = curve_fit(Wcostheta_func, phi_d_2H_Wcostheta_costheta_center[index[i]:index[i+1]], phi_d_2H_Wcostheta_results[index[i]:index[i+1]], p0=[0.0])
#     curve_fit_residuals = phi_d_2H_Wcostheta_results[index[i]:index[i+1]] - Wcostheta_func(phi_d_2H_Wcostheta_costheta_center[index[i]:index[i+1]], curve_fit_params[0])
#     reduced_chi2 = np.sum((curve_fit_residuals/phi_d_2H_Wcostheta_results_statserr[index[i]:index[i+1]])**2)/(len(phi_d_2H_Wcostheta_results[index[i]:index[i+1]])-1)
#     axs[i].plot(np.linspace(-1, 1, 100), Wcostheta_func(np.linspace(-1, 1, 100), curve_fit_params[0]), 'b--', label='Fit')
#     axs[i].text(-0.9, 0.95, r'$\rho^0_{00}=%.2f\pm%.2f$' % (curve_fit_params[0], np.sqrt(curve_fit_cov[0])), fontsize=10, color='b', ha='left', va='top')
#     axs[i].text(-0.9, 0.85, r'$\chi^2/\mathrm{dof}=%.2f$' % reduced_chi2, fontsize=10, color='b', ha='left', va='top')
#     axs[i].plot(np.linspace(-1, 1, 100), 0.75*(1-np.linspace(-1, 1, 100)**2), 'r--', label='SCHC+NPE')
# axs[0].set_xlim(-1, 1)
# axs[0].set_ylim(0, 1)
# axs[0].set_xticks(np.arange(-0.75, 0.9, 0.25))
# axs[0].set_ylabel(r'$W(\cos\vartheta)$')
# fig.suptitle(r"$d(\gamma, \phi d')$ normalized distribution of $\cos\vartheta$$")
# fig.supxlabel(r'$\cos\vartheta$')
# plt.legend()
# file_pdf.savefig()
# plt.close()

# #======================================================================phi_d_2H_Wphi======================================================================

# # Plot the results
# fig = plt.figure(figsize=(6*(len(index) - 1), 6), dpi=300)
# gs = fig.add_gridspec(1, len(index) - 1, wspace=0)
# axs = gs.subplots(sharex=True, sharey=True)
# for i in range(len(index) - 1):
#     axs[i].errorbar(phi_d_2H_Wphi_phi_center[index[i]:index[i+1]], phi_d_2H_Wphi_results[index[i]:index[i+1]], xerr=phi_d_2H_Wphi_phi_width[index[i]:index[i+1]], yerr=phi_d_2H_Wphi_results_statserr[index[i]:index[i+1]], fmt='k.', label='This work')
#     axs[i].set_title(r'$%.1f<\mathrm{E}_{\gamma}<%.1f\ \mathrm{GeV},\  %.1f<-t<%.1f\ \mathrm{GeV}^2$' % (phi_d_2H_Wphi_energy_low[index[i]], phi_d_2H_Wphi_energy_high[index[i]], phi_d_2H_Wphi_minust_low[index[i]], phi_d_2H_Wphi_minust_high[index[i]]))
#     curve_fit_params, curve_fit_cov = curve_fit(Wphi_func, phi_d_2H_Wphi_phi_center[index[i]:index[i+1]], phi_d_2H_Wphi_results[index[i]:index[i+1]], p0=[0.0])
#     curve_fit_residuals = phi_d_2H_Wphi_results[index[i]:index[i+1]] - Wphi_func(phi_d_2H_Wphi_phi_center[index[i]:index[i+1]], curve_fit_params[0])
#     reduced_chi2 = np.sum((curve_fit_residuals/phi_d_2H_Wphi_results_statserr[index[i]:index[i+1]])**2)/(len(phi_d_2H_Wphi_results[index[i]:index[i+1]])-1)
#     axs[i].plot(np.linspace(-180, 180, 360), Wphi_func(np.linspace(-180, 180, 360), curve_fit_params[0]), 'b--', label='Fit')
#     axs[i].text(-175, 1.95, r'$\rho^0_{1-1}=%.2f\pm%.2f$' % (curve_fit_params[0], np.sqrt(curve_fit_cov[0])), fontsize=10, color='b', ha='left', va='top')
#     axs[i].text(-175, 1.85, r'$\chi^2/\mathrm{dof}=%.2f$' % reduced_chi2, fontsize=10, color='b', ha='left', va='top')
#     axs[i].plot(np.linspace(-180, 180, 360), np.ones(360), 'r--', label='SCHC+NPE')
# axs[0].set_xlim(-180, 180)
# axs[0].set_ylim(0, 2)
# axs[0].set_xticks(np.arange(-120, 180, 60))
# axs[0].set_ylabel(r'$2\pi W(\varphi)$')
# fig.suptitle(r"$d(\gamma, \phi d')$ normalized distribution of $\varphi$")
# fig.supxlabel(r'$\varphi\ [\mathrm{deg}]$')
# plt.legend()
# file_pdf.savefig()
# plt.close()

# #======================================================================phi_d_2H_WPhi======================================================================

# # Read the bin edges
# phi_d_2H_WPhi_energy_low            = np.loadtxt('configs/bins_phi_d_WPhi.txt')[:,0]
# phi_d_2H_WPhi_energy_high           = np.loadtxt('configs/bins_phi_d_WPhi.txt')[:,1]
# phi_d_2H_WPhi_minust_low            = np.loadtxt('configs/bins_phi_d_WPhi.txt')[:,2]
# phi_d_2H_WPhi_minust_high           = np.loadtxt('configs/bins_phi_d_WPhi.txt')[:,3]
# phi_d_2H_WPhi_Phi_low               = np.loadtxt('configs/bins_phi_d_WPhi.txt')[:,4]
# phi_d_2H_WPhi_Phi_high              = np.loadtxt('configs/bins_phi_d_WPhi.txt')[:,5]
# phi_d_2H_WPhi_Phi_center            = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_WPhi.txt')[:,8]
# phi_d_2H_WPhi_Phi_width             = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_WPhi.txt')[:,9]

# # Read the yield numbers
# phi_d_2H_WPhi_yield_data            = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_WPhi.txt')[:,12]
# phi_d_2H_WPhi_yield_sim             = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_WPhi.txt')[:,12]
# phi_d_2H_WPhi_yield_tagged          = np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_WPhi.txt')[:,12]
# phi_d_2H_WPhi_yield_data_statserr   = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_WPhi.txt')[:,13]/phi_d_2H_WPhi_yield_data
# phi_d_2H_WPhi_yield_sim_statserr    = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_WPhi.txt')[:,13]/phi_d_2H_WPhi_yield_sim
# phi_d_2H_WPhi_yield_tagged_statserr = np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_WPhi.txt')[:,13]/phi_d_2H_WPhi_yield_tagged

# # Calculate the efficiency
# phi_d_2H_WPhi_efficiency            = phi_d_2H_WPhi_yield_sim/phi_d_2H_WPhi_yield_tagged
# phi_d_2H_WPhi_efficiency_statserr   = phi_d_2H_WPhi_efficiency*np.sqrt(phi_d_2H_WPhi_yield_sim_statserr**2 + phi_d_2H_WPhi_yield_tagged_statserr**2)

# # Calculate the results
# phi_d_2H_WPhi_results               = phi_d_2H_WPhi_yield_data/phi_d_2H_WPhi_efficiency  # raw results
# phi_d_2H_WPhi_results               = normalize_distribution(phi_d_2H_WPhi_results, phi_d_2H_WPhi_energy_low, phi_d_2H_WPhi_minust_low) # normalize to have the sum equal to 1
# phi_d_2H_WPhi_results               = 2*np.pi*phi_d_2H_WPhi_results/((phi_d_2H_WPhi_Phi_high - phi_d_2H_WPhi_Phi_low)/180*np.pi)  # normalize to have the integral equal to 2pi
# phi_d_2H_WPhi_results_statserr      = phi_d_2H_WPhi_results*np.sqrt(phi_d_2H_WPhi_yield_data_statserr**2 + phi_d_2H_WPhi_efficiency_statserr**2)

# # Find the indices for the different energy and t bins
# index = []
# for i in range(len(phi_d_2H_WPhi_results)):
#     if (i == 0):
#         index.append(i)
#     elif (i == len(phi_d_2H_WPhi_results) - 1):
#         index.append(i+1)
#     else:
#         if (phi_d_2H_WPhi_energy_low[i] != phi_d_2H_WPhi_energy_low[i-1]) or (phi_d_2H_WPhi_minust_low[i] != phi_d_2H_WPhi_minust_low[i-1]):
#             index.append(i)

# # Plot the data yield
# fig = plt.figure(figsize=(6*(len(index) - 1), 6), dpi=300)
# gs = fig.add_gridspec(1, len(index) - 1, wspace=0)
# axs = gs.subplots(sharex=True, sharey=True)
# for i in range(len(index) - 1):
#     axs[i].errorbar(phi_d_2H_WPhi_Phi_center[index[i]:index[i+1]], phi_d_2H_WPhi_yield_data[index[i]:index[i+1]], xerr=phi_d_2H_WPhi_Phi_width[index[i]:index[i+1]], yerr=phi_d_2H_WPhi_yield_data_statserr[index[i]:index[i+1]], fmt='k.', label='This work')
#     axs[i].set_title(r'$%.1f<\mathrm{E}_{\gamma}<%.1f\ \mathrm{GeV},\  %.1f<-t<%.1f\ \mathrm{GeV}^2$' % (phi_d_2H_WPhi_energy_low[index[i]], phi_d_2H_WPhi_energy_high[index[i]], phi_d_2H_WPhi_minust_low[index[i]], phi_d_2H_WPhi_minust_high[index[i]]))
# axs[0].set_xlim(-180, 180)
# axs[0].set_ylim(0, 150)
# axs[0].set_xticks(np.arange(-120, 180, 60))
# axs[0].set_ylabel(r'$\mathrm{Yield}$')
# fig.suptitle(r"$d(\gamma, \phi d')$ yield vs $\Phi$")
# fig.supxlabel(r'$\Phi\ [\mathrm{deg}]$')
# plt.legend()
# file_pdf.savefig()
# plt.close()

# # Plot the efficiency
# fig = plt.figure(figsize=(6*(len(index) - 1), 6), dpi=300)
# gs = fig.add_gridspec(1, len(index) - 1, wspace=0)
# axs = gs.subplots(sharex=True, sharey=True)
# for i in range(len(index) - 1):
#     axs[i].errorbar(phi_d_2H_WPhi_Phi_center[index[i]:index[i+1]], phi_d_2H_WPhi_efficiency[index[i]:index[i+1]], xerr=phi_d_2H_WPhi_Phi_width[index[i]:index[i+1]], yerr=phi_d_2H_WPhi_efficiency_statserr[index[i]:index[i+1]], fmt='k.', label='This work')
#     axs[i].set_title(r'$%.1f<\mathrm{E}_{\gamma}<%.1f\ \mathrm{GeV},\  %.1f<-t<%.1f\ \mathrm{GeV}^2$' % (phi_d_2H_WPhi_energy_low[index[i]], phi_d_2H_WPhi_energy_high[index[i]], phi_d_2H_WPhi_minust_low[index[i]], phi_d_2H_WPhi_minust_high[index[i]]))
# axs[0].set_xlim(-180, 180)
# axs[0].set_ylim(0, 0.1)
# axs[0].set_xticks(np.arange(-120, 180, 60))
# axs[0].set_ylabel(r'$\mathrm{Efficiency}$')
# fig.suptitle(r"$d(\gamma, \phi d')$ efficiency vs $\Phi$")
# fig.supxlabel(r'$\Phi\ [\mathrm{deg}]$')
# plt.legend()
# file_pdf.savefig()
# plt.close()

# # Plot the results
# fig = plt.figure(figsize=(6*(len(index) - 1), 6), dpi=300)
# gs = fig.add_gridspec(1, len(index) - 1, wspace=0)
# axs = gs.subplots(sharex=True, sharey=True)
# for i in range(len(index) - 1):
#     axs[i].errorbar(phi_d_2H_WPhi_Phi_center[index[i]:index[i+1]], phi_d_2H_WPhi_results[index[i]:index[i+1]], xerr=phi_d_2H_WPhi_Phi_width[index[i]:index[i+1]], yerr=phi_d_2H_WPhi_results_statserr[index[i]:index[i+1]], fmt='k.', label='This work')
#     axs[i].set_title(r'$%.1f<\mathrm{E}_{\gamma}<%.1f\ \mathrm{GeV},\  %.1f<-t<%.1f\ \mathrm{GeV}^2$' % (phi_d_2H_WPhi_energy_low[index[i]], phi_d_2H_WPhi_energy_high[index[i]], phi_d_2H_WPhi_minust_low[index[i]], phi_d_2H_WPhi_minust_high[index[i]]))
#     curve_fit_params, curve_fit_cov = curve_fit(Wphi_func, phi_d_2H_WPhi_Phi_center[index[i]:index[i+1]], phi_d_2H_WPhi_results[index[i]:index[i+1]], p0=[0.0])
#     curve_fit_residuals = phi_d_2H_WPhi_results[index[i]:index[i+1]] - Wphi_func(phi_d_2H_WPhi_Phi_center[index[i]:index[i+1]], curve_fit_params[0])
#     reduced_chi2 = np.sum((curve_fit_residuals/phi_d_2H_WPhi_results_statserr[index[i]:index[i+1]])**2)/(len(phi_d_2H_WPhi_results[index[i]:index[i+1]])-1)
#     axs[i].plot(np.linspace(-180, 180, 360), Wphi_func(np.linspace(-180, 180, 360), curve_fit_params[0]), 'b--', label='Fit')
#     axs[i].text(-175, 1.95, r'$\rho^1_{00}=%.2f\pm%.2f$' % (2*curve_fit_params[0]/3/0.3, np.sqrt(curve_fit_cov[0])), fontsize=10, color='b', ha='left', va='top')
#     axs[i].text(-175, 1.85, r'$\chi^2/\mathrm{dof}=%.2f$' % reduced_chi2, fontsize=10, color='b', ha='left', va='top')
#     axs[i].plot(np.linspace(-180, 180, 360), np.ones(360), 'r--', label='SCHC+NPE')
# axs[0].set_xlim(-180, 180)
# axs[0].set_ylim(0, 2)
# axs[0].set_xticks(np.arange(-120, 180, 60))
# axs[0].set_ylabel(r'$2\pi W(\Phi)$')
# fig.suptitle(r"$d(\gamma, \phi d')$ normalized distribution of $\Phi$")
# fig.supxlabel(r'$\Phi\ [\mathrm{deg}]$')
# plt.legend()
# file_pdf.savefig()
# plt.close()

# #======================================================================phi_d_2H_Wpsi======================================================================

# # Read the bin edges
# phi_d_2H_Wpsi_energy_low            = np.loadtxt('configs/bins_phi_d_Wpsi.txt')[:,0]
# phi_d_2H_Wpsi_energy_high           = np.loadtxt('configs/bins_phi_d_Wpsi.txt')[:,1]
# phi_d_2H_Wpsi_minust_low            = np.loadtxt('configs/bins_phi_d_Wpsi.txt')[:,2]
# phi_d_2H_Wpsi_minust_high           = np.loadtxt('configs/bins_phi_d_Wpsi.txt')[:,3]
# phi_d_2H_Wpsi_psi_low               = np.loadtxt('configs/bins_phi_d_Wpsi.txt')[:,4]
# phi_d_2H_Wpsi_psi_high              = np.loadtxt('configs/bins_phi_d_Wpsi.txt')[:,5]
# phi_d_2H_Wpsi_psi_center            = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_Wpsi.txt')[:,8]
# phi_d_2H_Wpsi_psi_width             = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_Wpsi.txt')[:,9]

# # Read the yield numbers
# phi_d_2H_Wpsi_yield_data            = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_Wpsi.txt')[:,12]
# phi_d_2H_Wpsi_yield_sim             = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_Wpsi.txt')[:,12]
# phi_d_2H_Wpsi_yield_tagged          = np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_Wpsi.txt')[:,12]
# phi_d_2H_Wpsi_yield_data_statserr   = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_Wpsi.txt')[:,13]/phi_d_2H_Wpsi_yield_data
# phi_d_2H_Wpsi_yield_sim_statserr    = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_Wpsi.txt')[:,13]/phi_d_2H_Wpsi_yield_sim
# phi_d_2H_Wpsi_yield_tagged_statserr = np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_Wpsi.txt')[:,13]/phi_d_2H_Wpsi_yield_tagged

# # Calculate the efficiency
# phi_d_2H_Wpsi_efficiency            = phi_d_2H_Wpsi_yield_sim/phi_d_2H_Wpsi_yield_tagged
# phi_d_2H_Wpsi_efficiency_statserr   = phi_d_2H_Wpsi_efficiency*np.sqrt(phi_d_2H_Wpsi_yield_sim_statserr**2 + phi_d_2H_Wpsi_yield_tagged_statserr**2)

# # Calculate the results
# phi_d_2H_Wpsi_results               = phi_d_2H_Wpsi_yield_data/phi_d_2H_Wpsi_efficiency  # raw results
# phi_d_2H_Wpsi_results               = normalize_distribution(phi_d_2H_Wpsi_results, phi_d_2H_Wpsi_energy_low, phi_d_2H_Wpsi_minust_low) # normalize to have the sum equal to 1
# phi_d_2H_Wpsi_results               = 2*np.pi*phi_d_2H_Wpsi_results/((phi_d_2H_Wpsi_psi_high - phi_d_2H_Wpsi_psi_low)/180*np.pi)  # normalize to have the integral equal to 2pi
# phi_d_2H_Wpsi_results_statserr      = phi_d_2H_Wpsi_results*np.sqrt(phi_d_2H_Wpsi_yield_data_statserr**2 + phi_d_2H_Wpsi_efficiency_statserr**2)

# # Find the indices for the different energy and t bins
# index = []
# for i in range(len(phi_d_2H_Wpsi_results)):
#     if (i == 0):
#         index.append(i)
#     elif (i == len(phi_d_2H_Wpsi_results) - 1):
#         index.append(i+1)
#     else:
#         if (phi_d_2H_Wpsi_energy_low[i] != phi_d_2H_Wpsi_energy_low[i-1]) or (phi_d_2H_Wpsi_minust_low[i] != phi_d_2H_Wpsi_minust_low[i-1]):
#             index.append(i)

# # Plot the data yield
# fig = plt.figure(figsize=(6*(len(index) - 1), 6), dpi=300)
# gs = fig.add_gridspec(1, len(index) - 1, wspace=0)
# axs = gs.subplots(sharex=True, sharey=True)
# for i in range(len(index) - 1):
#     axs[i].errorbar(phi_d_2H_Wpsi_psi_center[index[i]:index[i+1]], phi_d_2H_Wpsi_yield_data[index[i]:index[i+1]], xerr=phi_d_2H_Wpsi_psi_width[index[i]:index[i+1]], yerr=phi_d_2H_Wpsi_yield_data_statserr[index[i]:index[i+1]], fmt='k.', label='This work')
#     axs[i].set_title(r'$%.1f<\mathrm{E}_{\gamma}<%.1f\ \mathrm{GeV},\  %.1f<-t<%.1f\ \mathrm{GeV}^2$' % (phi_d_2H_Wpsi_energy_low[index[i]], phi_d_2H_Wpsi_energy_high[index[i]], phi_d_2H_Wpsi_minust_low[index[i]], phi_d_2H_Wpsi_minust_high[index[i]]))
# axs[0].set_xlim(-180, 180)
# axs[0].set_ylim(0, 150)
# axs[0].set_xticks(np.arange(-120, 180, 60))
# axs[0].set_ylabel(r'$\mathrm{Yield}$')
# fig.suptitle(r"$d(\gamma, \phi d')$ yield vs $\psi$")
# fig.supxlabel(r'$\psi\ [\mathrm{deg}]$')
# plt.legend()
# file_pdf.savefig()
# plt.close()

# # Plot the efficiency
# fig = plt.figure(figsize=(6*(len(index) - 1), 6), dpi=300)
# gs = fig.add_gridspec(1, len(index) - 1, wspace=0)
# axs = gs.subplots(sharex=True, sharey=True)
# for i in range(len(index) - 1):
#     axs[i].errorbar(phi_d_2H_Wpsi_psi_center[index[i]:index[i+1]], phi_d_2H_Wpsi_efficiency[index[i]:index[i+1]], xerr=phi_d_2H_Wpsi_psi_width[index[i]:index[i+1]], yerr=phi_d_2H_Wpsi_efficiency_statserr[index[i]:index[i+1]], fmt='k.', label='This work')
#     axs[i].set_title(r'$%.1f<\mathrm{E}_{\gamma}<%.1f\ \mathrm{GeV},\  %.1f<-t<%.1f\ \mathrm{GeV}^2$' % (phi_d_2H_Wpsi_energy_low[index[i]], phi_d_2H_Wpsi_energy_high[index[i]], phi_d_2H_Wpsi_minust_low[index[i]], phi_d_2H_Wpsi_minust_high[index[i]]))
# axs[0].set_xlim(-180, 180)
# axs[0].set_ylim(0, 0.1)
# axs[0].set_xticks(np.arange(-120, 180, 60))
# axs[0].set_ylabel(r'$\mathrm{Efficiency}$')
# fig.suptitle(r"$d(\gamma, \phi d')$ efficiency vs $\psi$")
# fig.supxlabel(r'$\psi\ [\mathrm{deg}]$')
# plt.legend()
# file_pdf.savefig()
# plt.close()

# # Plot the results
# fig = plt.figure(figsize=(6*(len(index) - 1), 6), dpi=300)
# gs = fig.add_gridspec(1, len(index) - 1, wspace=0)
# axs = gs.subplots(sharex=True, sharey=True)
# for i in range(len(index) - 1):
#     axs[i].errorbar(phi_d_2H_Wpsi_psi_center[index[i]:index[i+1]], phi_d_2H_Wpsi_results[index[i]:index[i+1]], xerr=phi_d_2H_Wpsi_psi_width[index[i]:index[i+1]], yerr=phi_d_2H_Wpsi_results_statserr[index[i]:index[i+1]], fmt='k.', label='This work')
#     axs[i].set_title(r'$%.1f<\mathrm{E}_{\gamma}<%.1f\ \mathrm{GeV},\  %.1f<-t<%.1f\ \mathrm{GeV}^2$' % (phi_d_2H_Wpsi_energy_low[index[i]], phi_d_2H_Wpsi_energy_high[index[i]], phi_d_2H_Wpsi_minust_low[index[i]], phi_d_2H_Wpsi_minust_high[index[i]]))
#     curve_fit_params, curve_fit_cov = curve_fit(Wphi_func, phi_d_2H_Wpsi_psi_center[index[i]:index[i+1]], phi_d_2H_Wpsi_results[index[i]:index[i+1]], p0=[0.0])
#     curve_fit_residuals = phi_d_2H_Wpsi_results[index[i]:index[i+1]] - Wphi_func(phi_d_2H_Wpsi_psi_center[index[i]:index[i+1]], curve_fit_params[0])
#     reduced_chi2 = np.sum((curve_fit_residuals/phi_d_2H_Wpsi_results_statserr[index[i]:index[i+1]])**2)/(len(phi_d_2H_Wpsi_results[index[i]:index[i+1]])-1)
#     axs[i].plot(np.linspace(-180, 180, 360), Wphi_func(np.linspace(-180, 180, 360), curve_fit_params[0]), 'b--', label='Fit')
#     axs[i].text(-175, 1.95, r'$\rho^1_{1-1}=%.2f\pm%.2f$' % (-curve_fit_params[0]/0.2, np.sqrt(curve_fit_cov[0])), fontsize=10, color='b', ha='left', va='top')
#     axs[i].text(-175, 1.85, r'$\chi^2/\mathrm{dof}=%.2f$' % reduced_chi2, fontsize=10, color='b', ha='left', va='top')
#     axs[i].plot(np.linspace(-180, 180, 360), np.ones(360)+0.2*np.cos(2*np.linspace(-180, 180, 360)/rad_to_deg), 'r--', label='SCHC+NPE')
# axs[0].set_xlim(-180, 180)
# axs[0].set_ylim(0, 2)
# axs[0].set_xticks(np.arange(-120, 180, 60))
# axs[0].set_ylabel(r'$2\pi W(\psi)$')
# fig.suptitle(r"$d(\gamma, \phi d')$ normalized distribution of $\psi$")
# fig.supxlabel(r'$\psi\ [\mathrm{deg}]$')
# plt.legend()
# file_pdf.savefig()
# plt.close()

# # Read the bin edges
# phi_d_2H_Wcostheta_energy_low            = np.loadtxt('configs/bins_phi_d_Wcostheta.txt')[:,0]
# phi_d_2H_Wcostheta_energy_high           = np.loadtxt('configs/bins_phi_d_Wcostheta.txt')[:,1]
# phi_d_2H_Wcostheta_minust_low            = np.loadtxt('configs/bins_phi_d_Wcostheta.txt')[:,2]
# phi_d_2H_Wcostheta_minust_high           = np.loadtxt('configs/bins_phi_d_Wcostheta.txt')[:,3]
# phi_d_2H_Wcostheta_costheta_low          = np.loadtxt('configs/bins_phi_d_Wcostheta.txt')[:,4]
# phi_d_2H_Wcostheta_costheta_high         = np.loadtxt('configs/bins_phi_d_Wcostheta.txt')[:,5]
# phi_d_2H_Wcostheta_costheta_center       = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wcostheta_sideband.txt')[:,8]
# phi_d_2H_Wcostheta_costheta_width        = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wcostheta_sideband.txt')[:,9]

# # Read the yield numbers
# phi_d_2H_Wcostheta_yield_data            = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wcostheta_sideband.txt')[:,12]
# phi_d_2H_Wcostheta_yield_data_statserr   = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wcostheta_sideband.txt')[:,13]

# # Calculate the results
# phi_d_2H_Wcostheta_results                  = phi_d_2H_Wcostheta_yield_data # raw results
# phi_d_2H_Wcostheta_results_statserr         = (1/np.sqrt(phi_d_2H_Wcostheta_yield_data_statserr**2) + 1/np.sqrt(np.sum(phi_d_2H_Wcostheta_yield_data_statserr**2)))*phi_d_2H_Wcostheta_results
# phi_d_2H_Wcostheta_results                  = normalize_distribution(phi_d_2H_Wcostheta_results, phi_d_2H_Wcostheta_results_statserr, phi_d_2H_Wcostheta_energy_low, phi_d_2H_Wcostheta_minust_low) # normalize to have the sum equal to 1
# phi_d_2H_Wcostheta_results                  = phi_d_2H_Wcostheta_results/(phi_d_2H_Wcostheta_costheta_high - phi_d_2H_Wcostheta_costheta_low)  # normalize to have the integral equal to 1

# # Find the indices for the different energy and t bins
# index = []
# for i in range(len(phi_d_2H_Wcostheta_results)):
#     if (i == 0):
#         index.append(i)
#     elif (i == len(phi_d_2H_Wcostheta_results) - 1):
#         index.append(i+1)
#     else:
#         if (phi_d_2H_Wcostheta_energy_low[i] != phi_d_2H_Wcostheta_energy_low[i-1]) or (phi_d_2H_Wcostheta_minust_low[i] != phi_d_2H_Wcostheta_minust_low[i-1]):
#             index.append(i)

# # Plot the results
# fig = plt.figure(figsize=(6*(len(index) - 1), 6), dpi=300)
# gs = fig.add_gridspec(1, len(index) - 1, wspace=0)
# axs = gs.subplots(sharex=True, sharey=True)
# for i in range(len(index) - 1):
#     axs[i].errorbar(phi_d_2H_Wcostheta_costheta_center[index[i]:index[i+1]], phi_d_2H_Wcostheta_results[index[i]:index[i+1]], xerr=phi_d_2H_Wcostheta_costheta_width[index[i]:index[i+1]], yerr=phi_d_2H_Wcostheta_results_statserr[index[i]:index[i+1]], fmt='k.', label='This work')
#     axs[i].set_title(r'$%.1f<\mathrm{E}_{\gamma}<%.1f\ \mathrm{GeV},\  %.1f<-t<%.1f\ \mathrm{GeV}^2$' % (phi_d_2H_Wcostheta_energy_low[index[i]], phi_d_2H_Wcostheta_energy_high[index[i]], phi_d_2H_Wcostheta_minust_low[index[i]], phi_d_2H_Wcostheta_minust_high[index[i]]))
#     axs[i].plot(np.linspace(-1, 1, 100), 0.75*(1-np.linspace(-1, 1, 100)**2), 'r--', label='SCHC+NPE')
# axs[0].set_xlim(-1, 1)
# axs[0].set_ylim(0, 1)
# axs[0].set_xticks(np.arange(-0.75, 0.9, 0.25))
# axs[0].set_ylabel(r'$W(\cos\vartheta)$')
# fig.suptitle(r"$d(\gamma, \phi d')$ normalized distribution of $\cos\vartheta$$")
# fig.supxlabel(r'$\cos\vartheta$')
# plt.legend()
# file_pdf.savefig()
# plt.close()

# # Read the bin edges
# phi_d_2H_Wpolphi_energy_low            = np.loadtxt('configs/bins_phi_d_Wpolphi.txt')[:,0]
# phi_d_2H_Wpolphi_energy_high           = np.loadtxt('configs/bins_phi_d_Wpolphi.txt')[:,1]
# phi_d_2H_Wpolphi_minust_low            = np.loadtxt('configs/bins_phi_d_Wpolphi.txt')[:,2]
# phi_d_2H_Wpolphi_minust_high           = np.loadtxt('configs/bins_phi_d_Wpolphi.txt')[:,3]
# phi_d_2H_Wpolphi_polphi_low               = np.loadtxt('configs/bins_phi_d_Wpolphi.txt')[:,4]
# phi_d_2H_Wpolphi_polphi_high              = np.loadtxt('configs/bins_phi_d_Wpolphi.txt')[:,5]
# phi_d_2H_Wpolphi_polphi_center            = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpolphi_sideband.txt')[:,8]
# phi_d_2H_Wpolphi_polphi_width             = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpolphi_sideband.txt')[:,9]

# # Read the yield numbers
# phi_d_2H_Wpolphi_yield_data            = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpolphi_sideband.txt')[:,12]
# phi_d_2H_Wpolphi_yield_data_statserr   = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpolphi_sideband.txt')[:,13]

# # Calculate the results
# phi_d_2H_Wpolphi_results               = phi_d_2H_Wpolphi_yield_data
# phi_d_2H_Wpolphi_results               = normalize_distribution(phi_d_2H_Wpolphi_results, phi_d_2H_Wpolphi_energy_low, phi_d_2H_Wpolphi_minust_low) # normalize to have the sum equal to 1
# phi_d_2H_Wpolphi_results               = 2*np.pi*phi_d_2H_Wpolphi_results/((phi_d_2H_Wpolphi_polphi_high - phi_d_2H_Wpolphi_polphi_low)/180*np.pi)  # normalize to have the integral equal to 2pi
# phi_d_2H_Wpolphi_results_statserr      = np.sqrt(1/phi_d_2H_Wpolphi_yield_data_statserr**2 + 1/np.sum(phi_d_2H_Wpolphi_yield_data_statserr**2))*phi_d_2H_Wpolphi_results

# # Find the indices for the different energy and t bins
# index = []
# for i in range(len(phi_d_2H_Wpolphi_results)):
#     if (i == 0):
#         index.append(i)
#     elif (i == len(phi_d_2H_Wpolphi_results) - 1):
#         index.append(i+1)
#     else:
#         if (phi_d_2H_Wpolphi_energy_low[i] != phi_d_2H_Wpolphi_energy_low[i-1]) or (phi_d_2H_Wpolphi_minust_low[i] != phi_d_2H_Wpolphi_minust_low[i-1]):
#             index.append(i)

# # Plot the results
# fig = plt.figure(figsize=(6*(len(index) - 1), 6), dpi=300)
# gs = fig.add_gridspec(1, len(index) - 1, wspace=0)
# axs = gs.subplots(sharex=True, sharey=True)
# for i in range(len(index) - 1):
#     axs[i].errorbar(phi_d_2H_Wpolphi_polphi_center[index[i]:index[i+1]], phi_d_2H_Wpolphi_results[index[i]:index[i+1]], xerr=phi_d_2H_Wpolphi_polphi_width[index[i]:index[i+1]], yerr=phi_d_2H_Wpolphi_results_statserr[index[i]:index[i+1]], fmt='k.', label='This work')
#     axs[i].set_title(r'$%.1f<\mathrm{E}_{\gamma}<%.1f\ \mathrm{GeV},\  %.1f<-t<%.1f\ \mathrm{GeV}^2$' % (phi_d_2H_Wpolphi_energy_low[index[i]], phi_d_2H_Wpolphi_energy_high[index[i]], phi_d_2H_Wpolphi_minust_low[index[i]], phi_d_2H_Wpolphi_minust_high[index[i]]))
# axs[0].set_xlim(-180, 180)
# axs[0].set_ylim(0, 2)
# axs[0].set_xticks(np.arange(-120, 180, 60))
# axs[0].set_ylabel(r'$2\pi W(\varphi)$')
# fig.suptitle(r"$d(\gamma, \phi d')$ normalized distribution of $\varphi$")
# fig.supxlabel(r'$\varphi\ [\mathrm{deg}]$')
# plt.legend()
# file_pdf.savefig()
# plt.close()

# # Read the bin edges
# phi_d_2H_Wdecayphi_energy_low            = np.loadtxt('configs/bins_phi_d_Wdecayphi.txt')[:,0]
# phi_d_2H_Wdecayphi_energy_high           = np.loadtxt('configs/bins_phi_d_Wdecayphi.txt')[:,1]
# phi_d_2H_Wdecayphi_minust_low            = np.loadtxt('configs/bins_phi_d_Wdecayphi.txt')[:,2]
# phi_d_2H_Wdecayphi_minust_high           = np.loadtxt('configs/bins_phi_d_Wdecayphi.txt')[:,3]
# phi_d_2H_Wdecayphi_decayphi_low               = np.loadtxt('configs/bins_phi_d_Wdecayphi.txt')[:,4]
# phi_d_2H_Wdecayphi_decayphi_high              = np.loadtxt('configs/bins_phi_d_Wdecayphi.txt')[:,5]
# phi_d_2H_Wdecayphi_decayphi_center            = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wdecayphi_sideband.txt')[:,8]
# phi_d_2H_Wdecayphi_decayphi_width             = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wdecayphi_sideband.txt')[:,9]

# # Read the yield numbers
# phi_d_2H_Wdecayphi_yield_data            = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wdecayphi_sideband.txt')[:,12]
# phi_d_2H_Wdecayphi_yield_data_statserr   = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wdecayphi_sideband.txt')[:,13]

# # Calculate the results
# phi_d_2H_Wdecayphi_results               = phi_d_2H_Wdecayphi_yield_data
# phi_d_2H_Wdecayphi_results               = normalize_distribution(phi_d_2H_Wdecayphi_results, phi_d_2H_Wdecayphi_energy_low, phi_d_2H_Wdecayphi_minust_low) # normalize to have the sum equal to 1
# phi_d_2H_Wdecayphi_results               = 2*np.pi*phi_d_2H_Wdecayphi_results/((phi_d_2H_Wdecayphi_decayphi_high - phi_d_2H_Wdecayphi_decayphi_low)/180*np.pi)  # normalize to have the integral equal to 2pi
# phi_d_2H_Wdecayphi_results_statserr      = (1/np.sqrt(phi_d_2H_Wdecayphi_yield_data_statserr**2) + 1/np.sqrt(np.sum(phi_d_2H_Wdecayphi_yield_data_statserr**2)))*phi_d_2H_Wdecayphi_results

# # Find the indices for the different energy and t bins
# index = []
# for i in range(len(phi_d_2H_Wdecayphi_results)):
#     if (i == 0):
#         index.append(i)
#     elif (i == len(phi_d_2H_Wdecayphi_results) - 1):
#         index.append(i+1)
#     else:
#         if (phi_d_2H_Wdecayphi_energy_low[i] != phi_d_2H_Wdecayphi_energy_low[i-1]) or (phi_d_2H_Wdecayphi_minust_low[i] != phi_d_2H_Wdecayphi_minust_low[i-1]):
#             index.append(i)

# # Plot the results
# fig = plt.figure(figsize=(6*(len(index) - 1), 6), dpi=300)
# gs = fig.add_gridspec(1, len(index) - 1, wspace=0)
# axs = gs.subplots(sharex=True, sharey=True)
# for i in range(len(index) - 1):
#     axs[i].errorbar(phi_d_2H_Wdecayphi_decayphi_center[index[i]:index[i+1]], phi_d_2H_Wdecayphi_results[index[i]:index[i+1]], xerr=phi_d_2H_Wdecayphi_decayphi_width[index[i]:index[i+1]], yerr=phi_d_2H_Wdecayphi_results_statserr[index[i]:index[i+1]], fmt='k.', label='This work')
#     axs[i].set_title(r'$%.1f<\mathrm{E}_{\gamma}<%.1f\ \mathrm{GeV},\  %.1f<-t<%.1f\ \mathrm{GeV}^2$' % (phi_d_2H_Wdecayphi_energy_low[index[i]], phi_d_2H_Wdecayphi_energy_high[index[i]], phi_d_2H_Wdecayphi_minust_low[index[i]], phi_d_2H_Wdecayphi_minust_high[index[i]]))
# axs[0].set_xlim(-180, 180)
# axs[0].set_ylim(0, 2)
# axs[0].set_xticks(np.arange(-120, 180, 60))
# axs[0].set_ylabel(r'$2\pi W(\varphi)$')
# fig.suptitle(r"$d(\gamma, \phi d')$ normalized distribution of $\varphi$")
# fig.supxlabel(r'$\varphi\ [\mathrm{deg}]$')
# plt.legend()
# file_pdf.savefig()
# plt.close()

# # Read the bin edges
# phi_d_2H_Wpsi_energy_low            = np.loadtxt('configs/bins_phi_d_Wpsi.txt')[:,0]
# phi_d_2H_Wpsi_energy_high           = np.loadtxt('configs/bins_phi_d_Wpsi.txt')[:,1]
# phi_d_2H_Wpsi_minust_low            = np.loadtxt('configs/bins_phi_d_Wpsi.txt')[:,2]
# phi_d_2H_Wpsi_minust_high           = np.loadtxt('configs/bins_phi_d_Wpsi.txt')[:,3]
# phi_d_2H_Wpsi_psi_low               = np.loadtxt('configs/bins_phi_d_Wpsi.txt')[:,4]
# phi_d_2H_Wpsi_psi_high              = np.loadtxt('configs/bins_phi_d_Wpsi.txt')[:,5]
# phi_d_2H_Wpsi_psi_center            = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpsi_sideband.txt')[:,8]
# phi_d_2H_Wpsi_psi_width             = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpsi_sideband.txt')[:,9]

# # Read the yield numbers
# phi_d_2H_Wpsi_yield_data            = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpsi_sideband.txt')[:,12]
# phi_d_2H_Wpsi_yield_data_statserr   = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpsi_sideband.txt')[:,13]

# # Calculate the results
# phi_d_2H_Wpsi_results               = phi_d_2H_Wpsi_yield_data
# phi_d_2H_Wpsi_results               = normalize_distribution(phi_d_2H_Wpsi_results, phi_d_2H_Wpsi_energy_low, phi_d_2H_Wpsi_minust_low) # normalize to have the sum equal to 1
# phi_d_2H_Wpsi_results               = 2*np.pi*phi_d_2H_Wpsi_results/((phi_d_2H_Wpsi_psi_high - phi_d_2H_Wpsi_psi_low)/180*np.pi)  # normalize to have the integral equal to 2pi
# phi_d_2H_Wpsi_results_statserr      = (1/np.sqrt(phi_d_2H_Wpsi_yield_data_statserr**2) + 1/np.sqrt(np.sum(phi_d_2H_Wpsi_yield_data_statserr**2)))*phi_d_2H_Wpsi_results

# # Find the indices for the different energy and t bins
# index = []
# for i in range(len(phi_d_2H_Wpsi_results)):
#     if (i == 0):
#         index.append(i)
#     elif (i == len(phi_d_2H_Wpsi_results) - 1):
#         index.append(i+1)
#     else:
#         if (phi_d_2H_Wpsi_energy_low[i] != phi_d_2H_Wpsi_energy_low[i-1]) or (phi_d_2H_Wpsi_minust_low[i] != phi_d_2H_Wpsi_minust_low[i-1]):
#             index.append(i)

# # Plot the results
# fig = plt.figure(figsize=(6*(len(index) - 1), 6), dpi=300)
# gs = fig.add_gridspec(1, len(index) - 1, wspace=0)
# axs = gs.subplots(sharex=True, sharey=True)
# for i in range(len(index) - 1):
#     axs[i].errorbar(phi_d_2H_Wpsi_psi_center[index[i]:index[i+1]], phi_d_2H_Wpsi_results[index[i]:index[i+1]], xerr=phi_d_2H_Wpsi_psi_width[index[i]:index[i+1]], yerr=phi_d_2H_Wpsi_results_statserr[index[i]:index[i+1]], fmt='k.', label='This work')
#     axs[i].set_title(r'$%.1f<\mathrm{E}_{\gamma}<%.1f\ \mathrm{GeV},\  %.1f<-t<%.1f\ \mathrm{GeV}^2$' % (phi_d_2H_Wpsi_energy_low[index[i]], phi_d_2H_Wpsi_energy_high[index[i]], phi_d_2H_Wpsi_minust_low[index[i]], phi_d_2H_Wpsi_minust_high[index[i]]))
# axs[0].set_xlim(-180, 180)
# axs[0].set_ylim(0, 2)
# axs[0].set_xticks(np.arange(-120, 180, 60))
# axs[0].set_ylabel(r'$2\pi W(\psi)$')
# fig.suptitle(r"$d(\gamma, \phi d')$ normalized distribution of $\psi$")
# fig.supxlabel(r'$\psi\ [\mathrm{deg}]$')
# plt.legend()
# file_pdf.savefig()
# plt.close()

file_pdf.close()