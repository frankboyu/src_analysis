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
                sigma_this                  = results_statserr[index:i+1]
                sigma_rest                  = np.sqrt(np.sum(results_statserr[index:i+1]**2) - sigma_this**2)
                results[index:i+1]          = results[index:i+1]/N_total
                results_statserr[index:i+1] = np.sqrt(N_rest**2 * sigma_this**2 + N_this**2 * sigma_rest**2)/N_total**2
        else:
            if (energy_bins[i] != energy_bins[i-1]) or (t_bins[i] != t_bins[i-1]):
                N_total                     = np.sum(results[index:i])
                N_this                      = results[index:i]
                N_rest                      = np.sum(results[index:i]) - N_this
                sigma_this                  = results_statserr[index:i]
                sigma_rest                  = np.sqrt(np.sum(results_statserr[index:i]**2) - sigma_this**2)
                results[index:i]            = results[index:i]/N_total
                results_statserr[index:i]   = np.sqrt(N_rest**2 * sigma_this**2 + N_this**2 * sigma_rest**2)/N_total**2
                index = i
    return results, results_statserr

def legend_without_duplicate_labels(ax):
    handles, labels = ax.get_legend_handles_labels()
    unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
    ax.legend(*zip(*unique))

################################################################# LISTS OF LABELS ###################################################################################

dsdt_label_list = [r'$E_\gamma$=5.8-7.8 GeV', r'$E_\gamma$=7.8-8.8 GeV', r'$E_\gamma$=8.8-10.8 GeV']
dsdt_color_list = ['r', 'g', 'b']
W_label_list    = [ r'$E_\gamma$=5.8-7.8 GeV, $-t<0.6 \rm{ GeV}^2$', \
                    r'$E_\gamma$=5.8-7.8 GeV, $-t>0.6 \rm{ GeV}^2$', \
                    r'$E_\gamma$=7.8-8.8 GeV, $-t<0.6 \rm{ GeV}^2$', \
                    r'$E_\gamma$=7.8-8.8 GeV, $-t>0.6 \rm{ GeV}^2$', \
                    r'$E_\gamma$=8.8-10.8 GeV, $-t<0.6 \rm{ GeV}^2$', \
                    r'$E_\gamma$=8.8-10.8 GeV, $-t>0.6 \rm{ GeV}^2$']
W_color_list    = ['r', 'orange', 'green', 'cyan', 'blue', 'purple']

xlabel_list = [r'$-t[GeV^2/c]$',    r'$\cos\vartheta$',     r'$\varphi$[deg]',  r'$\Phi$[deg]', r'$\psi$[deg]']
title_list  = [r'$d\sigma /dt$',    r'$W(\cos\vartheta)$',  r'$W(\varphi$)',    r'$W(\Phi$)',   r'$W(\psi$)']
xmin_list   = [0,                   -1,                     -180,               -180,           -180]
xmax_list   = [2,                   1,                      180,                180,            180]

color_options   = (['r', 'orange', 'black', 'cyan', 'blue', 'green', 'yellow', 'purple', 'brown', 'pink', 'gray', 'olive', 'navy'])
cut_name_list   = (['dEdx cut', 'Missing p- cut', 'Chi2/NDF cut', 'Momentum cut', 'Theta cut', 'Vertex Z cut', 'Vertex R cut', 'Beam accid subtraction', 'Combo accid subtraction', 'Fit max', 'Fit width', 'Fit background model', 'Fit signal model', 'Simulation parameter a1', 'Simulation parameter b1', 'Simulation parameter a2', 'Simulation parameter b2'])
cut_class_list  = ([         1,                1,              1,              1,           1,              1,              1,                        2,                         2,         3,           3,                      3,                  3,                         4,                         4,                         4,                         4])
p2p_list        = (['Total p2p', 'Event selection', 'Accidental subtraction', 'Fitting procedure', 'Simulation model'])

nominal_list        = (['nominal'])
dedx_list           = (['dEdx_1.50',                'dEdx_1.75',                'dEdx_2.50',                'dEdx_3.00'])
pminus_list         = (['misspminus_0.0150',        'misspminus_0.0175',        'misspminus_0.0250',        'misspminus_0.0300'])
chisquared_list     = (['chisquared_4.50',          'chisquared_4.75',          'chisquared_5.50',          'chisquared_6.00'])
momentum_list       = (['momentum_0.350',           'momentum_0.375',           'momentum_0.425',           'momentum_0.450'])
theta_list          = (['theta_1.90',               'theta_1.95',               'theta_2.05',               'theta_2.10'])
vertexZ_list        = (['vertexZ_13.50',            'vertexZ_13.75',            'vertexZ_14.25',            'vertexZ_14.50'])
vertexR_list        = (['vertexR_0.50',             'vertexR_0.75',             'vertexR_1.25',             'vertexR_1.50'])
beamaccid_list      = (['beamaccid_3',              'beamaccid_5',              'beamaccid_4out'])
comboaccid_list     = (['comboaccid_all',           'comboaccid_none'])
fitmax_list         = (['fitmax_1.06',              'fitmax_1.07',              'fitmax_1.09',              'fitmax_1.10'])
fitwidth_list       = (['fitwidth_0.0040',          'fitwidth_0.0048',          'fitwidth_0.0060',          'fitwidth_0.0075'])
fitbkg_list         = (['fitbkg_fulllinear',        'fitbkg_quadratic',         'fitbkg_fullquadratic',     'fitbkg_phenomenological'])
fitsig_list         = (['fitsig_noBL',              'fitsig_nonrel'])
simweight_a1_list   = (['simweight_syst_a1_-1.0',   'simweight_syst_a1_-0.5',   'simweight_syst_a1_0.5',    'simweight_syst_a1_1.0'])
simweight_b1_list   = (['simweight_syst_b1_-1.0',   'simweight_syst_b1_-0.5',   'simweight_syst_b1_0.5',    'simweight_syst_b1_1.0'])
simweight_a2_list   = (['simweight_syst_a2_-1.0',   'simweight_syst_a2_-0.5',   'simweight_syst_a2_0.5',    'simweight_syst_a2_1.0'])
simweight_b2_list   = (['simweight_syst_b2_-1.0',   'simweight_syst_b2_-0.5',   'simweight_syst_b2_0.5',    'simweight_syst_b2_1.0'])

variation_name_list = []
cut_index_list      = []
color_list          = []
for i, this_list in enumerate([nominal_list, dedx_list, pminus_list, chisquared_list, momentum_list, theta_list, vertexZ_list, vertexR_list, beamaccid_list, comboaccid_list, fitmax_list, fitwidth_list, fitbkg_list, fitsig_list, simweight_a1_list, simweight_b1_list, simweight_a2_list, simweight_b2_list]):
    for j, this_label in enumerate(this_list):
            variation_name_list.append(this_label)
            cut_index_list.append(i-1)
            color_list.append(color_options[j])

################################################################# READ THE NUMBERS ###################################################################################

phi_d_2H_dsdt_energy_center                 = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_dsdt_nominal.txt')[:,0]
phi_d_2H_dsdt_energy_width                  = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_dsdt_nominal.txt')[:,1]
phi_d_2H_dsdt_energy_low                    = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_dsdt_nominal.txt')[:,2]
phi_d_2H_dsdt_energy_high                   = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_dsdt_nominal.txt')[:,3]
phi_d_2H_dsdt_minust_center                 = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_dsdt_nominal.txt')[:,4]
phi_d_2H_dsdt_minust_width                  = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_dsdt_nominal.txt')[:,5]
phi_d_2H_dsdt_minust_low                    = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_dsdt_nominal.txt')[:,6]
phi_d_2H_dsdt_minust_high                   = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_dsdt_nominal.txt')[:,7]

phi_d_2H_Wcostheta_energy_center            = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wcostheta_nominal.txt')[:,0]
phi_d_2H_Wcostheta_energy_width             = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wcostheta_nominal.txt')[:,1]
phi_d_2H_Wcostheta_energy_low               = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wcostheta_nominal.txt')[:,2]
phi_d_2H_Wcostheta_energy_high              = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wcostheta_nominal.txt')[:,3]
phi_d_2H_Wcostheta_minust_center            = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wcostheta_nominal.txt')[:,4]
phi_d_2H_Wcostheta_minust_width             = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wcostheta_nominal.txt')[:,5]
phi_d_2H_Wcostheta_minust_low               = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wcostheta_nominal.txt')[:,6]
phi_d_2H_Wcostheta_minust_high              = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wcostheta_nominal.txt')[:,7]
phi_d_2H_Wcostheta_costheta_center          = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wcostheta_nominal.txt')[:,8]
phi_d_2H_Wcostheta_costheta_width           = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wcostheta_nominal.txt')[:,9]
phi_d_2H_Wcostheta_costheta_low             = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wcostheta_nominal.txt')[:,10]
phi_d_2H_Wcostheta_costheta_high            = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wcostheta_nominal.txt')[:,11]

phi_d_2H_Wdecayphi_energy_center            = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wdecayphi_nominal.txt')[:,0]
phi_d_2H_Wdecayphi_energy_width             = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wdecayphi_nominal.txt')[:,1]
phi_d_2H_Wdecayphi_energy_low               = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wdecayphi_nominal.txt')[:,2]
phi_d_2H_Wdecayphi_energy_high              = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wdecayphi_nominal.txt')[:,3]
phi_d_2H_Wdecayphi_minust_center            = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wdecayphi_nominal.txt')[:,4]
phi_d_2H_Wdecayphi_minust_width             = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wdecayphi_nominal.txt')[:,5]
phi_d_2H_Wdecayphi_minust_low               = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wdecayphi_nominal.txt')[:,6]
phi_d_2H_Wdecayphi_minust_high              = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wdecayphi_nominal.txt')[:,7]
phi_d_2H_Wdecayphi_decayphi_center          = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wdecayphi_nominal.txt')[:,8]
phi_d_2H_Wdecayphi_decayphi_width           = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wdecayphi_nominal.txt')[:,9]
phi_d_2H_Wdecayphi_decayphi_low             = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wdecayphi_nominal.txt')[:,10]
phi_d_2H_Wdecayphi_decayphi_high            = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wdecayphi_nominal.txt')[:,11]

phi_d_2H_Wpolphi_energy_center              = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpolphi_nominal.txt')[:,0]
phi_d_2H_Wpolphi_energy_width               = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpolphi_nominal.txt')[:,1]
phi_d_2H_Wpolphi_energy_low                 = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpolphi_nominal.txt')[:,2]
phi_d_2H_Wpolphi_energy_high                = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpolphi_nominal.txt')[:,3]
phi_d_2H_Wpolphi_minust_center              = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpolphi_nominal.txt')[:,4]
phi_d_2H_Wpolphi_minust_width               = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpolphi_nominal.txt')[:,5]
phi_d_2H_Wpolphi_minust_low                 = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpolphi_nominal.txt')[:,6]
phi_d_2H_Wpolphi_minust_high                = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpolphi_nominal.txt')[:,7]
phi_d_2H_Wpolphi_polphi_center              = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpolphi_nominal.txt')[:,8]
phi_d_2H_Wpolphi_polphi_width               = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpolphi_nominal.txt')[:,9]
phi_d_2H_Wpolphi_polphi_low                 = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpolphi_nominal.txt')[:,10]
phi_d_2H_Wpolphi_polphi_high                = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpolphi_nominal.txt')[:,11]

phi_d_2H_Wpsi_energy_center                 = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpsi_nominal.txt')[:,0]
phi_d_2H_Wpsi_energy_width                  = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpsi_nominal.txt')[:,1]
phi_d_2H_Wpsi_energy_low                    = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpsi_nominal.txt')[:,2]
phi_d_2H_Wpsi_energy_high                   = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpsi_nominal.txt')[:,3]
phi_d_2H_Wpsi_minust_center                 = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpsi_nominal.txt')[:,4]
phi_d_2H_Wpsi_minust_width                  = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpsi_nominal.txt')[:,5]
phi_d_2H_Wpsi_minust_low                    = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpsi_nominal.txt')[:,6]
phi_d_2H_Wpsi_minust_high                   = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpsi_nominal.txt')[:,7]
phi_d_2H_Wpsi_psi_center                    = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpsi_nominal.txt')[:,8]
phi_d_2H_Wpsi_psi_width                     = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpsi_nominal.txt')[:,9]
phi_d_2H_Wpsi_psi_low                       = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpsi_nominal.txt')[:,10]
phi_d_2H_Wpsi_psi_high                      = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpsi_nominal.txt')[:,11]

# Read the yield numbers
phi_d_2H_dsdt_yield_data                    = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_dsdt_nominal.txt')[:,8]
phi_d_2H_dsdt_yield_data_statserr           = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_dsdt_nominal.txt')[:,9]
phi_d_2H_dsdt_yield_sim                     = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_flat_dsdt_nominal.txt')[:,8]
phi_d_2H_dsdt_yield_sim_statserr            = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_flat_dsdt_nominal.txt')[:,9]
phi_d_2H_dsdt_yield_tagged                  = np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_flat_dsdt_nominal.txt')[:,8]
phi_d_2H_dsdt_yield_tagged_statserr         = np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_flat_dsdt_nominal.txt')[:,9]
phi_d_2H_dsdt_efficiency                    = phi_d_2H_dsdt_yield_sim/phi_d_2H_dsdt_yield_tagged
phi_d_2H_dsdt_efficiency_statserr           = phi_d_2H_dsdt_efficiency*np.sqrt((phi_d_2H_dsdt_yield_sim_statserr/phi_d_2H_dsdt_yield_sim)**2 + (phi_d_2H_dsdt_yield_tagged_statserr/phi_d_2H_dsdt_yield_tagged)**2)
phi_d_2H_dsdt_results                       = phi_d_2H_dsdt_yield_data/phi_d_2H_dsdt_efficiency/lumi(phi_d_2H_dsdt_energy_low, phi_d_2H_dsdt_energy_high, 28)/(phi_d_2H_dsdt_minust_high-phi_d_2H_dsdt_minust_low)/0.489/1000
phi_d_2H_dsdt_results_statserr              = phi_d_2H_dsdt_results*np.sqrt((phi_d_2H_dsdt_yield_data_statserr/phi_d_2H_dsdt_yield_data)**2 + (phi_d_2H_dsdt_efficiency_statserr/phi_d_2H_dsdt_efficiency)**2)

phi_d_2H_Wcostheta_yield_data               = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wcostheta_nominal.txt')[:,12]
phi_d_2H_Wcostheta_yield_data_statserr      = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wcostheta_nominal.txt')[:,13]
phi_d_2H_Wcostheta_yield_sim                = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_flat_Wcostheta_nominal.txt')[:,12]
phi_d_2H_Wcostheta_yield_sim_statserr       = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_flat_Wcostheta_nominal.txt')[:,13]
phi_d_2H_Wcostheta_yield_tagged             = np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_flat_Wcostheta_nominal.txt')[:,12]
phi_d_2H_Wcostheta_yield_tagged_statserr    = np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_flat_Wcostheta_nominal.txt')[:,13]
phi_d_2H_Wcostheta_efficiency               = phi_d_2H_Wcostheta_yield_sim/phi_d_2H_Wcostheta_yield_tagged
phi_d_2H_Wcostheta_efficiency_statserr      = phi_d_2H_Wcostheta_efficiency*np.sqrt((phi_d_2H_Wcostheta_yield_sim_statserr/phi_d_2H_Wcostheta_yield_sim)**2 + (phi_d_2H_Wcostheta_yield_tagged_statserr/phi_d_2H_Wcostheta_yield_tagged)**2)
phi_d_2H_Wcostheta_results                  = phi_d_2H_Wcostheta_yield_data/phi_d_2H_Wcostheta_efficiency  # raw results
phi_d_2H_Wcostheta_results_statserr         = phi_d_2H_Wcostheta_results*np.sqrt((phi_d_2H_Wcostheta_yield_data_statserr/phi_d_2H_Wcostheta_yield_data)**2 + (phi_d_2H_Wcostheta_efficiency_statserr/phi_d_2H_Wcostheta_efficiency)**2)
phi_d_2H_Wcostheta_results, \
phi_d_2H_Wcostheta_results_statserr         = normalize_distribution(phi_d_2H_Wcostheta_results, phi_d_2H_Wcostheta_results_statserr, phi_d_2H_Wcostheta_energy_low, phi_d_2H_Wcostheta_minust_low) # normalize to have the sum equal to 1
phi_d_2H_Wcostheta_results                  = phi_d_2H_Wcostheta_results/(phi_d_2H_Wcostheta_costheta_high - phi_d_2H_Wcostheta_costheta_low)  # normalize to have the integral equal to 1
phi_d_2H_Wcostheta_results_statserr         = phi_d_2H_Wcostheta_results_statserr/(phi_d_2H_Wcostheta_costheta_high - phi_d_2H_Wcostheta_costheta_low)  # set the proper statistical uncertainties

phi_d_2H_Wdecayphi_yield_data               = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wdecayphi_nominal.txt')[:,12]
phi_d_2H_Wdecayphi_yield_data_statserr      = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wdecayphi_nominal.txt')[:,13]
phi_d_2H_Wdecayphi_yield_sim                = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_flat_Wdecayphi_nominal.txt')[:,12]
phi_d_2H_Wdecayphi_yield_sim_statserr       = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_flat_Wdecayphi_nominal.txt')[:,13]
phi_d_2H_Wdecayphi_yield_tagged             = np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_flat_Wdecayphi_nominal.txt')[:,12]
phi_d_2H_Wdecayphi_yield_tagged_statserr    = np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_flat_Wdecayphi_nominal.txt')[:,13]
phi_d_2H_Wdecayphi_efficiency               = phi_d_2H_Wdecayphi_yield_sim/phi_d_2H_Wdecayphi_yield_tagged
phi_d_2H_Wdecayphi_efficiency_statserr      = phi_d_2H_Wdecayphi_efficiency*np.sqrt((phi_d_2H_Wdecayphi_yield_sim_statserr/phi_d_2H_Wdecayphi_yield_sim)**2 + (phi_d_2H_Wdecayphi_yield_tagged_statserr/phi_d_2H_Wdecayphi_yield_tagged)**2)
phi_d_2H_Wdecayphi_results                  = phi_d_2H_Wdecayphi_yield_data/phi_d_2H_Wdecayphi_efficiency  # raw results
phi_d_2H_Wdecayphi_results_statserr         = phi_d_2H_Wdecayphi_results*np.sqrt((phi_d_2H_Wdecayphi_yield_data_statserr/phi_d_2H_Wdecayphi_yield_data)**2 + (phi_d_2H_Wdecayphi_efficiency_statserr/phi_d_2H_Wdecayphi_efficiency)**2)
phi_d_2H_Wdecayphi_results, \
phi_d_2H_Wdecayphi_results_statserr         = normalize_distribution(phi_d_2H_Wdecayphi_results, phi_d_2H_Wdecayphi_results_statserr, phi_d_2H_Wdecayphi_energy_low, phi_d_2H_Wdecayphi_minust_low) # normalize to have the sum equal to 1
phi_d_2H_Wdecayphi_results                  = 2*np.pi*phi_d_2H_Wdecayphi_results/((phi_d_2H_Wdecayphi_decayphi_high - phi_d_2H_Wdecayphi_decayphi_low)/180*np.pi)  # normalize to have the integral equal to 2pi
phi_d_2H_Wdecayphi_results_statserr         = 2*np.pi*phi_d_2H_Wdecayphi_results_statserr/((phi_d_2H_Wdecayphi_decayphi_high - phi_d_2H_Wdecayphi_decayphi_low)/180*np.pi)  # set the proper statistical uncertainties

phi_d_2H_Wpolphi_yield_data                 = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpolphi_nominal.txt')[:,12]
phi_d_2H_Wpolphi_yield_data_statserr        = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpolphi_nominal.txt')[:,13]
phi_d_2H_Wpolphi_yield_sim                  = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_flat_Wpolphi_nominal.txt')[:,12]
phi_d_2H_Wpolphi_yield_sim_statserr         = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_flat_Wpolphi_nominal.txt')[:,13]
phi_d_2H_Wpolphi_yield_tagged               = np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_flat_Wpolphi_nominal.txt')[:,12]
phi_d_2H_Wpolphi_yield_tagged_statserr      = np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_flat_Wpolphi_nominal.txt')[:,13]
phi_d_2H_Wpolphi_efficiency                 = phi_d_2H_Wpolphi_yield_sim/phi_d_2H_Wpolphi_yield_tagged
phi_d_2H_Wpolphi_efficiency_statserr        = phi_d_2H_Wpolphi_efficiency*np.sqrt((phi_d_2H_Wpolphi_yield_sim_statserr/phi_d_2H_Wpolphi_yield_sim)**2 + (phi_d_2H_Wpolphi_yield_tagged_statserr/phi_d_2H_Wpolphi_yield_tagged)**2)
phi_d_2H_Wpolphi_results                    = phi_d_2H_Wpolphi_yield_data/phi_d_2H_Wpolphi_efficiency  # raw results
phi_d_2H_Wpolphi_results_statserr           = phi_d_2H_Wpolphi_results*np.sqrt((phi_d_2H_Wpolphi_yield_data_statserr/phi_d_2H_Wpolphi_yield_data)**2 + (phi_d_2H_Wpolphi_efficiency_statserr/phi_d_2H_Wpolphi_efficiency)**2)
phi_d_2H_Wpolphi_results, \
phi_d_2H_Wpolphi_results_statserr           = normalize_distribution(phi_d_2H_Wpolphi_results, phi_d_2H_Wpolphi_results_statserr, phi_d_2H_Wpolphi_energy_low, phi_d_2H_Wpolphi_minust_low) # normalize to have the sum equal to 1
phi_d_2H_Wpolphi_results                    = 2*np.pi*phi_d_2H_Wpolphi_results/((phi_d_2H_Wpolphi_polphi_high - phi_d_2H_Wpolphi_polphi_low)/180*np.pi)  # normalize to have the integral equal to 2pi
phi_d_2H_Wpolphi_results_statserr           = 2*np.pi*phi_d_2H_Wpolphi_results_statserr/((phi_d_2H_Wpolphi_polphi_high - phi_d_2H_Wpolphi_polphi_low)/180*np.pi)  # set the proper statistical uncertainties

phi_d_2H_Wpsi_yield_data                    = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpsi_nominal.txt')[:,12]
phi_d_2H_Wpsi_yield_data_statserr           = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpsi_nominal.txt')[:,13]
phi_d_2H_Wpsi_yield_sim                     = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_flat_Wpsi_nominal.txt')[:,12]
phi_d_2H_Wpsi_yield_sim_statserr            = np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_flat_Wpsi_nominal.txt')[:,13]
phi_d_2H_Wpsi_yield_tagged                  = np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_flat_Wpsi_nominal.txt')[:,12]
phi_d_2H_Wpsi_yield_tagged_statserr         = np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_flat_Wpsi_nominal.txt')[:,13]
phi_d_2H_Wpsi_efficiency                    = phi_d_2H_Wpsi_yield_sim/phi_d_2H_Wpsi_yield_tagged
phi_d_2H_Wpsi_efficiency_statserr           = phi_d_2H_Wpsi_efficiency*np.sqrt((phi_d_2H_Wpsi_yield_sim_statserr/phi_d_2H_Wpsi_yield_sim)**2 + (phi_d_2H_Wpsi_yield_tagged_statserr/phi_d_2H_Wpsi_yield_tagged)**2)
phi_d_2H_Wpsi_results                       = phi_d_2H_Wpsi_yield_data/phi_d_2H_Wpsi_efficiency  # raw results
phi_d_2H_Wpsi_results_statserr              = phi_d_2H_Wpsi_results*np.sqrt((phi_d_2H_Wpsi_yield_data_statserr/phi_d_2H_Wpsi_yield_data)**2 + (phi_d_2H_Wpsi_efficiency_statserr/phi_d_2H_Wpsi_efficiency)**2)
phi_d_2H_Wpsi_results, \
phi_d_2H_Wpsi_results_statserr              = normalize_distribution(phi_d_2H_Wpsi_results, phi_d_2H_Wpsi_results_statserr, phi_d_2H_Wpsi_energy_low, phi_d_2H_Wpsi_minust_low) # normalize to have the sum equal to 1
phi_d_2H_Wpsi_results                       = 2*np.pi*phi_d_2H_Wpsi_results/((phi_d_2H_Wpsi_psi_high - phi_d_2H_Wpsi_psi_low)/180*np.pi)  # normalize to have the integral equal to 2pi
phi_d_2H_Wpsi_results_statserr              = 2*np.pi*phi_d_2H_Wpsi_results_statserr/((phi_d_2H_Wpsi_psi_high - phi_d_2H_Wpsi_psi_low)/180*np.pi)  # set the proper statistical uncertainties

for tag in variation_name_list:
    if (tag == 'nominal'):
        continue

    if (tag.find('simweight') != -1 or tag == 'fitsig_relBWsim'):
        phi_d_2H_dsdt_yield_data                    = np.vstack((phi_d_2H_dsdt_yield_data,                      np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_dsdt_nominal.txt')[:,8]))
        phi_d_2H_dsdt_yield_data_statserr           = np.vstack((phi_d_2H_dsdt_yield_data_statserr,             np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_dsdt_nominal.txt')[:,9]))
        phi_d_2H_Wcostheta_yield_data               = np.vstack((phi_d_2H_Wcostheta_yield_data,                 np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wcostheta_nominal.txt')[:,12]))
        phi_d_2H_Wcostheta_yield_data_statserr      = np.vstack((phi_d_2H_Wcostheta_yield_data_statserr,        np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wcostheta_nominal.txt')[:,13]))
        phi_d_2H_Wdecayphi_yield_data               = np.vstack((phi_d_2H_Wdecayphi_yield_data,                 np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wdecayphi_nominal.txt')[:,12]))
        phi_d_2H_Wdecayphi_yield_data_statserr      = np.vstack((phi_d_2H_Wdecayphi_yield_data_statserr,        np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wdecayphi_nominal.txt')[:,13]))
        phi_d_2H_Wpolphi_yield_data                 = np.vstack((phi_d_2H_Wpolphi_yield_data,                   np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpolphi_nominal.txt')[:,12]))
        phi_d_2H_Wpolphi_yield_data_statserr        = np.vstack((phi_d_2H_Wpolphi_yield_data_statserr,          np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpolphi_nominal.txt')[:,13]))
        phi_d_2H_Wpsi_yield_data                    = np.vstack((phi_d_2H_Wpsi_yield_data,                      np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpsi_nominal.txt')[:,12]))
        phi_d_2H_Wpsi_yield_data_statserr           = np.vstack((phi_d_2H_Wpsi_yield_data_statserr,             np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpsi_nominal.txt')[:,13]))
    else:
        phi_d_2H_dsdt_yield_data                    = np.vstack((phi_d_2H_dsdt_yield_data,                      np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_dsdt_'+tag+'.txt')[:,8]))
        phi_d_2H_dsdt_yield_data_statserr           = np.vstack((phi_d_2H_dsdt_yield_data_statserr,             np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_dsdt_'+tag+'.txt')[:,9]))
        phi_d_2H_Wcostheta_yield_data               = np.vstack((phi_d_2H_Wcostheta_yield_data,                 np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wcostheta_'+tag+'.txt')[:,12]))
        phi_d_2H_Wcostheta_yield_data_statserr      = np.vstack((phi_d_2H_Wcostheta_yield_data_statserr,        np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wcostheta_'+tag+'.txt')[:,13]))
        phi_d_2H_Wdecayphi_yield_data               = np.vstack((phi_d_2H_Wdecayphi_yield_data,                 np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wdecayphi_'+tag+'.txt')[:,12]))
        phi_d_2H_Wdecayphi_yield_data_statserr      = np.vstack((phi_d_2H_Wdecayphi_yield_data_statserr,        np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wdecayphi_'+tag+'.txt')[:,13]))
        phi_d_2H_Wpolphi_yield_data                 = np.vstack((phi_d_2H_Wpolphi_yield_data,                   np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpolphi_'+tag+'.txt')[:,12]))
        phi_d_2H_Wpolphi_yield_data_statserr        = np.vstack((phi_d_2H_Wpolphi_yield_data_statserr,          np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpolphi_'+tag+'.txt')[:,13]))
        phi_d_2H_Wpsi_yield_data                    = np.vstack((phi_d_2H_Wpsi_yield_data,                      np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpsi_'+tag+'.txt')[:,12]))
        phi_d_2H_Wpsi_yield_data_statserr           = np.vstack((phi_d_2H_Wpsi_yield_data_statserr,             np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_data_2H_ver12_Wpsi_'+tag+'.txt')[:,13]))

    if (tag.find('fit') != -1 and tag != 'fitsig_relBWsim'):
        phi_d_2H_dsdt_yield_sim                     = np.vstack((phi_d_2H_dsdt_yield_sim,                       np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_flat_dsdt_nominal.txt')[:,8]))
        phi_d_2H_dsdt_yield_sim_statserr            = np.vstack((phi_d_2H_dsdt_yield_sim_statserr,              np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_flat_dsdt_nominal.txt')[:,9]))
        phi_d_2H_Wcostheta_yield_sim                = np.vstack((phi_d_2H_Wcostheta_yield_sim,                  np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_flat_Wcostheta_nominal.txt')[:,12]))
        phi_d_2H_Wcostheta_yield_sim_statserr       = np.vstack((phi_d_2H_Wcostheta_yield_sim_statserr,         np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_flat_Wcostheta_nominal.txt')[:,13]))
        phi_d_2H_Wdecayphi_yield_sim                = np.vstack((phi_d_2H_Wdecayphi_yield_sim,                  np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_flat_Wdecayphi_nominal.txt')[:,12]))
        phi_d_2H_Wdecayphi_yield_sim_statserr       = np.vstack((phi_d_2H_Wdecayphi_yield_sim_statserr,         np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_flat_Wdecayphi_nominal.txt')[:,13]))
        phi_d_2H_Wpolphi_yield_sim                  = np.vstack((phi_d_2H_Wpolphi_yield_sim,                    np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_flat_Wpolphi_nominal.txt')[:,12]))
        phi_d_2H_Wpolphi_yield_sim_statserr         = np.vstack((phi_d_2H_Wpolphi_yield_sim_statserr,           np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_flat_Wpolphi_nominal.txt')[:,13]))
        phi_d_2H_Wpsi_yield_sim                     = np.vstack((phi_d_2H_Wpsi_yield_sim,                       np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_flat_Wpsi_nominal.txt')[:,12]))
        phi_d_2H_Wpsi_yield_sim_statserr            = np.vstack((phi_d_2H_Wpsi_yield_sim_statserr,              np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_flat_Wpsi_nominal.txt')[:,13]))
    else:
        phi_d_2H_dsdt_yield_sim                     = np.vstack((phi_d_2H_dsdt_yield_sim,                       np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_flat_dsdt_'+tag+'.txt')[:,8]))
        phi_d_2H_dsdt_yield_sim_statserr            = np.vstack((phi_d_2H_dsdt_yield_sim_statserr,              np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_flat_dsdt_'+tag+'.txt')[:,9]))
        phi_d_2H_Wcostheta_yield_sim                = np.vstack((phi_d_2H_Wcostheta_yield_sim,                  np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_flat_Wcostheta_'+tag+'.txt')[:,12]))
        phi_d_2H_Wcostheta_yield_sim_statserr       = np.vstack((phi_d_2H_Wcostheta_yield_sim_statserr,         np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_flat_Wcostheta_'+tag+'.txt')[:,13]))
        phi_d_2H_Wdecayphi_yield_sim                = np.vstack((phi_d_2H_Wdecayphi_yield_sim,                  np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_flat_Wdecayphi_'+tag+'.txt')[:,12]))
        phi_d_2H_Wdecayphi_yield_sim_statserr       = np.vstack((phi_d_2H_Wdecayphi_yield_sim_statserr,         np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_flat_Wdecayphi_'+tag+'.txt')[:,13]))
        phi_d_2H_Wpolphi_yield_sim                  = np.vstack((phi_d_2H_Wpolphi_yield_sim,                    np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_flat_Wpolphi_'+tag+'.txt')[:,12]))
        phi_d_2H_Wpolphi_yield_sim_statserr         = np.vstack((phi_d_2H_Wpolphi_yield_sim_statserr,           np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_flat_Wpolphi_'+tag+'.txt')[:,13]))
        phi_d_2H_Wpsi_yield_sim                     = np.vstack((phi_d_2H_Wpsi_yield_sim,                       np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_flat_Wpsi_'+tag+'.txt')[:,12]))
        phi_d_2H_Wpsi_yield_sim_statserr            = np.vstack((phi_d_2H_Wpsi_yield_sim_statserr,              np.loadtxt('output/yield_phi_d/yield_phi_d_recon_exc_sim_2H_ver12_flat_Wpsi_'+tag+'.txt')[:,13]))

    if (tag.find('simweight') != -1):
        phi_d_2H_dsdt_yield_tagged                  = np.vstack((phi_d_2H_dsdt_yield_tagged,                    np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_flat_dsdt_'+tag+'.txt')[:,8]))
        phi_d_2H_dsdt_yield_tagged_statserr         = np.vstack((phi_d_2H_dsdt_yield_tagged_statserr,           np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_flat_dsdt_'+tag+'.txt')[:,9]))
        phi_d_2H_Wcostheta_yield_tagged             = np.vstack((phi_d_2H_Wcostheta_yield_tagged,               np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_flat_Wcostheta_nominal.txt')[:,12]))
        phi_d_2H_Wcostheta_yield_tagged_statserr    = np.vstack((phi_d_2H_Wcostheta_yield_tagged_statserr,      np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_flat_Wcostheta_nominal.txt')[:,13]))
        phi_d_2H_Wdecayphi_yield_tagged             = np.vstack((phi_d_2H_Wdecayphi_yield_tagged,               np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_flat_Wdecayphi_nominal.txt')[:,12]))
        phi_d_2H_Wdecayphi_yield_tagged_statserr    = np.vstack((phi_d_2H_Wdecayphi_yield_tagged_statserr,      np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_flat_Wdecayphi_nominal.txt')[:,13]))
        phi_d_2H_Wpolphi_yield_tagged               = np.vstack((phi_d_2H_Wpolphi_yield_tagged,                 np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_flat_Wpolphi_nominal.txt')[:,12]))
        phi_d_2H_Wpolphi_yield_tagged_statserr      = np.vstack((phi_d_2H_Wpolphi_yield_tagged_statserr,        np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_flat_Wpolphi_nominal.txt')[:,13]))
        phi_d_2H_Wpsi_yield_tagged                  = np.vstack((phi_d_2H_Wpsi_yield_tagged,                    np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_flat_Wpsi_nominal.txt')[:,12]))
        phi_d_2H_Wpsi_yield_tagged_statserr         = np.vstack((phi_d_2H_Wpsi_yield_tagged_statserr,           np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_flat_Wpsi_nominal.txt')[:,13]))
    else:
        phi_d_2H_dsdt_yield_tagged                  = np.vstack((phi_d_2H_dsdt_yield_tagged,                    np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_flat_dsdt_nominal.txt')[:,8]))
        phi_d_2H_dsdt_yield_tagged_statserr         = np.vstack((phi_d_2H_dsdt_yield_tagged_statserr,           np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_flat_dsdt_nominal.txt')[:,9]))
        phi_d_2H_Wcostheta_yield_tagged             = np.vstack((phi_d_2H_Wcostheta_yield_tagged,               np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_flat_Wcostheta_nominal.txt')[:,12]))
        phi_d_2H_Wcostheta_yield_tagged_statserr    = np.vstack((phi_d_2H_Wcostheta_yield_tagged_statserr,      np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_flat_Wcostheta_nominal.txt')[:,13]))
        phi_d_2H_Wdecayphi_yield_tagged             = np.vstack((phi_d_2H_Wdecayphi_yield_tagged,               np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_flat_Wdecayphi_nominal.txt')[:,12]))
        phi_d_2H_Wdecayphi_yield_tagged_statserr    = np.vstack((phi_d_2H_Wdecayphi_yield_tagged_statserr,      np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_flat_Wdecayphi_nominal.txt')[:,13]))
        phi_d_2H_Wpolphi_yield_tagged               = np.vstack((phi_d_2H_Wpolphi_yield_tagged,                 np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_flat_Wpolphi_nominal.txt')[:,12]))
        phi_d_2H_Wpolphi_yield_tagged_statserr      = np.vstack((phi_d_2H_Wpolphi_yield_tagged_statserr,        np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_flat_Wpolphi_nominal.txt')[:,13]))
        phi_d_2H_Wpsi_yield_tagged                  = np.vstack((phi_d_2H_Wpsi_yield_tagged,                    np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_flat_Wpsi_nominal.txt')[:,12]))
        phi_d_2H_Wpsi_yield_tagged_statserr         = np.vstack((phi_d_2H_Wpsi_yield_tagged_statserr,           np.loadtxt('output/yield_phi_d/yield_phi_d_thrown_exc_tagged_2H_ver12_flat_Wpsi_nominal.txt')[:,13]))

    phi_d_2H_dsdt_efficiency                        = np.vstack((phi_d_2H_dsdt_efficiency,                      phi_d_2H_dsdt_yield_sim[-1,:]/phi_d_2H_dsdt_yield_tagged[-1,:]))
    phi_d_2H_dsdt_efficiency_statserr               = np.vstack((phi_d_2H_dsdt_efficiency_statserr,             phi_d_2H_dsdt_efficiency[-1,:]*np.sqrt((phi_d_2H_dsdt_yield_sim_statserr[-1,:]/phi_d_2H_dsdt_yield_sim[-1,:])**2 + (phi_d_2H_dsdt_yield_tagged_statserr[-1,:]/phi_d_2H_dsdt_yield_tagged[-1,:])**2)))
    phi_d_2H_Wcostheta_efficiency                   = np.vstack((phi_d_2H_Wcostheta_efficiency,                 phi_d_2H_Wcostheta_yield_sim[-1,:]/phi_d_2H_Wcostheta_yield_tagged[-1,:]))
    phi_d_2H_Wcostheta_efficiency_statserr          = np.vstack((phi_d_2H_Wcostheta_efficiency_statserr,        phi_d_2H_Wcostheta_efficiency[-1,:]*np.sqrt((phi_d_2H_Wcostheta_yield_sim_statserr[-1,:]/phi_d_2H_Wcostheta_yield_sim[-1,:])**2 + (phi_d_2H_Wcostheta_yield_tagged_statserr[-1,:]/phi_d_2H_Wcostheta_yield_tagged[-1,:])**2)))
    phi_d_2H_Wdecayphi_efficiency                   = np.vstack((phi_d_2H_Wdecayphi_efficiency,                 phi_d_2H_Wdecayphi_yield_sim[-1,:]/phi_d_2H_Wdecayphi_yield_tagged[-1,:]))
    phi_d_2H_Wdecayphi_efficiency_statserr          = np.vstack((phi_d_2H_Wdecayphi_efficiency_statserr,        phi_d_2H_Wdecayphi_efficiency[-1,:]*np.sqrt((phi_d_2H_Wdecayphi_yield_sim_statserr[-1,:]/phi_d_2H_Wdecayphi_yield_sim[-1,:])**2 + (phi_d_2H_Wdecayphi_yield_tagged_statserr[-1,:]/phi_d_2H_Wdecayphi_yield_tagged[-1,:])**2)))
    phi_d_2H_Wpolphi_efficiency                     = np.vstack((phi_d_2H_Wpolphi_efficiency,                   phi_d_2H_Wpolphi_yield_sim[-1,:]/phi_d_2H_Wpolphi_yield_tagged[-1,:]))
    phi_d_2H_Wpolphi_efficiency_statserr            = np.vstack((phi_d_2H_Wpolphi_efficiency_statserr,          phi_d_2H_Wpolphi_efficiency[-1,:]*np.sqrt((phi_d_2H_Wpolphi_yield_sim_statserr[-1,:]/phi_d_2H_Wpolphi_yield_sim[-1,:])**2 + (phi_d_2H_Wpolphi_yield_tagged_statserr[-1,:]/phi_d_2H_Wpolphi_yield_tagged[-1,:])**2)))
    phi_d_2H_Wpsi_efficiency                        = np.vstack((phi_d_2H_Wpsi_efficiency,                      phi_d_2H_Wpsi_yield_sim[-1,:]/phi_d_2H_Wpsi_yield_tagged[-1,:]))
    phi_d_2H_Wpsi_efficiency_statserr               = np.vstack((phi_d_2H_Wpsi_efficiency_statserr,             phi_d_2H_Wpsi_efficiency[-1,:]*np.sqrt((phi_d_2H_Wpsi_yield_sim_statserr[-1,:]/phi_d_2H_Wpsi_yield_sim[-1,:])**2 + (phi_d_2H_Wpsi_yield_tagged_statserr[-1,:]/phi_d_2H_Wpsi_yield_tagged[-1,:])**2)))

    phi_d_2H_dsdt_results                           = np.vstack((phi_d_2H_dsdt_results,                         phi_d_2H_dsdt_yield_data[-1,:]/phi_d_2H_dsdt_efficiency[-1,:]/lumi(phi_d_2H_dsdt_energy_low, phi_d_2H_dsdt_energy_high, 28)/(phi_d_2H_dsdt_minust_high-phi_d_2H_dsdt_minust_low)/0.489/1000))
    phi_d_2H_dsdt_results_statserr                  = np.vstack((phi_d_2H_dsdt_results_statserr,                phi_d_2H_dsdt_results[-1,:]*np.sqrt((phi_d_2H_dsdt_yield_data_statserr[-1,:]/phi_d_2H_dsdt_yield_data[-1,:])**2 + (phi_d_2H_dsdt_efficiency_statserr[-1,:]/phi_d_2H_dsdt_efficiency[-1,:])**2)))
    phi_d_2H_Wcostheta_results                      = np.vstack((phi_d_2H_Wcostheta_results,                    phi_d_2H_Wcostheta_yield_data[-1,:]/phi_d_2H_Wcostheta_efficiency[-1,:]))
    phi_d_2H_Wcostheta_results_statserr             = np.vstack((phi_d_2H_Wcostheta_results_statserr,           phi_d_2H_Wcostheta_results[-1,:]*np.sqrt((phi_d_2H_Wcostheta_yield_data_statserr[-1,:]/phi_d_2H_Wcostheta_yield_data[-1,:])**2 + (phi_d_2H_Wcostheta_efficiency_statserr[-1,:]/phi_d_2H_Wcostheta_efficiency[-1,:])**2)))
    phi_d_2H_Wcostheta_results[-1,:], \
    phi_d_2H_Wcostheta_results_statserr[-1,:]       = normalize_distribution(phi_d_2H_Wcostheta_results[-1,:],  phi_d_2H_Wcostheta_results_statserr[-1,:], phi_d_2H_Wcostheta_energy_low, phi_d_2H_Wcostheta_minust_low)
    phi_d_2H_Wcostheta_results[-1,:]                = phi_d_2H_Wcostheta_results[-1,:]/                         (phi_d_2H_Wcostheta_costheta_high - phi_d_2H_Wcostheta_costheta_low)
    phi_d_2H_Wcostheta_results_statserr[-1,:]       = phi_d_2H_Wcostheta_results_statserr[-1,:]/                (phi_d_2H_Wcostheta_costheta_high - phi_d_2H_Wcostheta_costheta_low)
    phi_d_2H_Wdecayphi_results                      = np.vstack((phi_d_2H_Wdecayphi_results,                    phi_d_2H_Wdecayphi_yield_data[-1,:]/phi_d_2H_Wdecayphi_efficiency[-1,:]))
    phi_d_2H_Wdecayphi_results_statserr             = np.vstack((phi_d_2H_Wdecayphi_results_statserr,           phi_d_2H_Wdecayphi_results[-1,:]*np.sqrt((phi_d_2H_Wdecayphi_yield_data_statserr[-1,:]/phi_d_2H_Wdecayphi_yield_data[-1,:])**2 + (phi_d_2H_Wdecayphi_efficiency_statserr[-1,:]/phi_d_2H_Wdecayphi_efficiency[-1,:])**2)))
    phi_d_2H_Wdecayphi_results[-1,:], \
    phi_d_2H_Wdecayphi_results_statserr[-1,:]       = normalize_distribution(phi_d_2H_Wdecayphi_results[-1,:],  phi_d_2H_Wdecayphi_results_statserr[-1,:], phi_d_2H_Wdecayphi_energy_low, phi_d_2H_Wdecayphi_minust_low)
    phi_d_2H_Wdecayphi_results[-1,:]                = 2*np.pi*phi_d_2H_Wdecayphi_results[-1,:]/                 ((phi_d_2H_Wdecayphi_decayphi_high - phi_d_2H_Wdecayphi_decayphi_low)/180*np.pi)
    phi_d_2H_Wdecayphi_results_statserr[-1,:]       = 2*np.pi*phi_d_2H_Wdecayphi_results_statserr[-1,:]/        ((phi_d_2H_Wdecayphi_decayphi_high - phi_d_2H_Wdecayphi_decayphi_low)/180*np.pi)
    phi_d_2H_Wpolphi_results                        = np.vstack((phi_d_2H_Wpolphi_results,                      phi_d_2H_Wpolphi_yield_data[-1,:]/phi_d_2H_Wpolphi_efficiency[-1,:]))
    phi_d_2H_Wpolphi_results_statserr               = np.vstack((phi_d_2H_Wpolphi_results_statserr,             phi_d_2H_Wpolphi_results[-1,:]*np.sqrt((phi_d_2H_Wpolphi_yield_data_statserr[-1,:]/phi_d_2H_Wpolphi_yield_data[-1,:])**2 + (phi_d_2H_Wpolphi_efficiency_statserr[-1,:]/phi_d_2H_Wpolphi_efficiency[-1,:])**2)))
    phi_d_2H_Wpolphi_results[-1,:], \
    phi_d_2H_Wpolphi_results_statserr[-1,:]         = normalize_distribution(phi_d_2H_Wpolphi_results[-1,:],    phi_d_2H_Wpolphi_results_statserr[-1,:], phi_d_2H_Wpolphi_energy_low, phi_d_2H_Wpolphi_minust_low)
    phi_d_2H_Wpolphi_results[-1,:]                  = 2*np.pi*phi_d_2H_Wpolphi_results[-1,:]/                   ((phi_d_2H_Wpolphi_polphi_high - phi_d_2H_Wpolphi_polphi_low)/180*np.pi)
    phi_d_2H_Wpolphi_results_statserr[-1,:]         = 2*np.pi*phi_d_2H_Wpolphi_results_statserr[-1,:]/          ((phi_d_2H_Wpolphi_polphi_high - phi_d_2H_Wpolphi_polphi_low)/180*np.pi)
    phi_d_2H_Wpsi_results                           = np.vstack((phi_d_2H_Wpsi_results,                         phi_d_2H_Wpsi_yield_data[-1,:]/phi_d_2H_Wpsi_efficiency[-1,:]))
    phi_d_2H_Wpsi_results_statserr                  = np.vstack((phi_d_2H_Wpsi_results_statserr,                phi_d_2H_Wpsi_results[-1,:]*np.sqrt((phi_d_2H_Wpsi_yield_data_statserr[-1,:]/phi_d_2H_Wpsi_yield_data[-1,:])**2 + (phi_d_2H_Wpsi_efficiency_statserr[-1,:]/phi_d_2H_Wpsi_efficiency[-1,:])**2)))
    phi_d_2H_Wpsi_results[-1,:], \
    phi_d_2H_Wpsi_results_statserr[-1,:]            = normalize_distribution(phi_d_2H_Wpsi_results[-1,:],       phi_d_2H_Wpsi_results_statserr[-1,:], phi_d_2H_Wpsi_energy_low, phi_d_2H_Wpsi_minust_low)
    phi_d_2H_Wpsi_results[-1,:]                     = 2*np.pi*phi_d_2H_Wpsi_results[-1,:]/                      ((phi_d_2H_Wpsi_psi_high - phi_d_2H_Wpsi_psi_low)/180*np.pi)
    phi_d_2H_Wpsi_results_statserr[-1,:]            = 2*np.pi*phi_d_2H_Wpsi_results_statserr[-1,:]/             ((phi_d_2H_Wpsi_psi_high - phi_d_2H_Wpsi_psi_low)/180*np.pi)

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

Wcostheta_index = []
for i in range(len(phi_d_2H_Wcostheta_yield_data[0,:])):
    if (i == 0):
        Wcostheta_index.append(i)
    elif (i == len(phi_d_2H_Wcostheta_yield_data[0,:]) - 1):
        Wcostheta_index.append(i+1)
    else:
        if ((phi_d_2H_Wcostheta_energy_low[i] != phi_d_2H_Wcostheta_energy_low[i-1]) or (phi_d_2H_Wcostheta_minust_low[i] != phi_d_2H_Wcostheta_minust_low[i-1])):
            Wcostheta_index.append(i)

Wdecayphi_index = []
for i in range(len(phi_d_2H_Wdecayphi_yield_data[0,:])):
    if (i == 0):
        Wdecayphi_index.append(i)
    elif (i == len(phi_d_2H_Wdecayphi_yield_data[0,:]) - 1):
        Wdecayphi_index.append(i+1)
    else:
        if ((phi_d_2H_Wdecayphi_energy_low[i] != phi_d_2H_Wdecayphi_energy_low[i-1]) or (phi_d_2H_Wdecayphi_minust_low[i] != phi_d_2H_Wdecayphi_minust_low[i-1])):
            Wdecayphi_index.append(i)

Wpolphi_index = []
for i in range(len(phi_d_2H_Wpolphi_yield_data[0,:])):
    if (i == 0):
        Wpolphi_index.append(i)
    elif (i == len(phi_d_2H_Wpolphi_yield_data[0,:]) - 1):
        Wpolphi_index.append(i+1)
    else:
        if ((phi_d_2H_Wpolphi_energy_low[i] != phi_d_2H_Wpolphi_energy_low[i-1]) or (phi_d_2H_Wpolphi_minust_low[i] != phi_d_2H_Wpolphi_minust_low[i-1])):
            Wpolphi_index.append(i)

Wpsi_index = []
for i in range(len(phi_d_2H_Wpsi_yield_data[0,:])):
    if (i == 0):
        Wpsi_index.append(i)
    elif (i == len(phi_d_2H_Wpsi_yield_data[0,:]) - 1):
        Wpsi_index.append(i+1)
    else:
        if ((phi_d_2H_Wpsi_energy_low[i] != phi_d_2H_Wpsi_energy_low[i-1]) or (phi_d_2H_Wpsi_minust_low[i] != phi_d_2H_Wpsi_minust_low[i-1])):
            Wpsi_index.append(i)

################################################################# DATA YIELD ##############################################################################################################################

print('Plotting data yield...')

fig = plt.figure(figsize=(16, 20), dpi=300)
gs = fig.add_gridspec(3, 4)
ax1 = fig.add_subplot(gs[0, 0:2])
ax2 = fig.add_subplot(gs[0, 2:4])
ax3 = fig.add_subplot(gs[1, 0:2])
ax4 = fig.add_subplot(gs[1, 2:4])
ax5 = fig.add_subplot(gs[2, 1:3])
axs = [ax1, ax2, ax3, ax4, ax5]

plt.subplots_adjust(wspace=0.5, hspace=0.2)

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
axs[0].set_ylim(0, 300)
axs[0].legend()

for i in range(len(Wcostheta_index)-1):
    axs[1].errorbar((phi_d_2H_Wcostheta_costheta_high[Wcostheta_index[i]:Wcostheta_index[i+1]] + phi_d_2H_Wcostheta_costheta_low[Wcostheta_index[i]:Wcostheta_index[i+1]])/2, \
                    phi_d_2H_Wcostheta_yield_data[0,Wcostheta_index[i]:Wcostheta_index[i+1]], \
                    xerr=(phi_d_2H_Wcostheta_costheta_high[Wcostheta_index[i]:Wcostheta_index[i+1]] - phi_d_2H_Wcostheta_costheta_low[Wcostheta_index[i]:Wcostheta_index[i+1]])/2, \
                    yerr=phi_d_2H_Wcostheta_yield_data_statserr[0,Wcostheta_index[i]:Wcostheta_index[i+1]], \
                    fmt='.', color=W_color_list[i], label=W_label_list[i])
axs[1].set_title(r'$W(\cos \vartheta)$')
axs[1].set_xlabel(r'$\cos\vartheta$')
axs[1].set_ylabel(r'$Y_{\rm data}$')
axs[1].set_xlim(-1, 1)
axs[1].set_ylim(0, 1000)
axs[1].legend(loc='upper left')

for i in range(len(Wdecayphi_index)-1):
    axs[2].errorbar((phi_d_2H_Wdecayphi_decayphi_high[Wdecayphi_index[i]:Wdecayphi_index[i+1]] + phi_d_2H_Wdecayphi_decayphi_low[Wdecayphi_index[i]:Wdecayphi_index[i+1]])/2, \
                    phi_d_2H_Wdecayphi_yield_data[0,Wdecayphi_index[i]:Wdecayphi_index[i+1]], \
                    xerr=(phi_d_2H_Wdecayphi_decayphi_high[Wdecayphi_index[i]:Wdecayphi_index[i+1]] - phi_d_2H_Wdecayphi_decayphi_low[Wdecayphi_index[i]:Wdecayphi_index[i+1]])/2, \
                    yerr=phi_d_2H_Wdecayphi_yield_data_statserr[0,Wdecayphi_index[i]:Wdecayphi_index[i+1]], \
                    fmt='.', color=W_color_list[i], label=W_label_list[i])
axs[2].set_title(r'$W(\varphi)$')
axs[2].set_xlabel(r'$\varphi$ (deg)')
axs[2].set_ylabel(r'$Y_{\rm data}$')
axs[2].set_xlim(-180, 180)
axs[2].set_ylim(0, 800)

for i in range(len(Wpolphi_index)-1):
    axs[3].errorbar((phi_d_2H_Wpolphi_polphi_high[Wpolphi_index[i]:Wpolphi_index[i+1]] + phi_d_2H_Wpolphi_polphi_low[Wpolphi_index[i]:Wpolphi_index[i+1]])/2, \
                    phi_d_2H_Wpolphi_yield_data[0,Wpolphi_index[i]:Wpolphi_index[i+1]], \
                    xerr=(phi_d_2H_Wpolphi_polphi_high[Wpolphi_index[i]:Wpolphi_index[i+1]] - phi_d_2H_Wpolphi_polphi_low[Wpolphi_index[i]:Wpolphi_index[i+1]])/2, \
                    yerr=phi_d_2H_Wpolphi_yield_data_statserr[0,Wpolphi_index[i]:Wpolphi_index[i+1]], \
                    fmt='.', color=W_color_list[i], label=W_label_list[i])
axs[3].set_title(r'$W(\Phi)$')
axs[3].set_xlabel(r'$\Phi$ (deg)')
axs[3].set_ylabel(r'$Y_{\rm data}$')
axs[3].set_xlim(-180, 180)
axs[3].set_ylim(0, 800)

for i in range(len(Wpsi_index)-1):
    axs[4].errorbar((phi_d_2H_Wpsi_psi_high[Wpsi_index[i]:Wpsi_index[i+1]] + phi_d_2H_Wpsi_psi_low[Wpsi_index[i]:Wpsi_index[i+1]])/2, \
                    phi_d_2H_Wpsi_yield_data[0,Wpsi_index[i]:Wpsi_index[i+1]], \
                    xerr=(phi_d_2H_Wpsi_psi_high[Wpsi_index[i]:Wpsi_index[i+1]] - phi_d_2H_Wpsi_psi_low[Wpsi_index[i]:Wpsi_index[i+1]])/2, \
                    yerr=phi_d_2H_Wpsi_yield_data_statserr[0,Wpsi_index[i]:Wpsi_index[i+1]], \
                    fmt='.', color=W_color_list[i], label=W_label_list[i])
axs[4].set_title(r'$W(\psi)$')
axs[4].set_xlabel(r'$\psi$ (deg)')
axs[4].set_ylabel(r'$Y_{\rm data}$')
axs[4].set_xlim(-180, 180)
axs[4].set_ylim(0, 800)

file_pdf.savefig()
plt.close()

################################################################# EFFICIENCY #############################################################################################################################

print('Plotting efficiency...')

fig = plt.figure(figsize=(16, 20), dpi=300)
gs = fig.add_gridspec(3, 4)
ax1 = fig.add_subplot(gs[0, 0:2])
ax2 = fig.add_subplot(gs[0, 2:4])
ax3 = fig.add_subplot(gs[1, 0:2])
ax4 = fig.add_subplot(gs[1, 2:4])
ax5 = fig.add_subplot(gs[2, 1:3])
axs = [ax1, ax2, ax3, ax4, ax5]

plt.subplots_adjust(wspace=0.5, hspace=0.2)

for i in range(len(dsdt_index)-1):
    axs[0].errorbar((phi_d_2H_dsdt_minust_high[dsdt_index[i]:dsdt_index[i+1]] + phi_d_2H_dsdt_minust_low[dsdt_index[i]:dsdt_index[i+1]])/2, \
                    phi_d_2H_dsdt_efficiency[0,dsdt_index[i]:dsdt_index[i+1]], \
                    xerr=(phi_d_2H_dsdt_minust_high[dsdt_index[i]:dsdt_index[i+1]] - phi_d_2H_dsdt_minust_low[dsdt_index[i]:dsdt_index[i+1]])/2, \
                    yerr=phi_d_2H_dsdt_efficiency_statserr[0,dsdt_index[i]:dsdt_index[i+1]], \
                    fmt=dsdt_color_list[i]+'.', label=dsdt_label_list[i])
axs[0].set_title(r'$d\sigma / dt$')
axs[0].set_xlabel(r'$-t[GeV^2/c]$')
axs[0].set_ylabel(r'$\epsilon$')
axs[0].set_xlim(0, 2)
axs[0].set_ylim(0, 0.5)
axs[0].legend()

for i in range(len(Wcostheta_index)-1):
    axs[1].errorbar((phi_d_2H_Wcostheta_costheta_high[Wcostheta_index[i]:Wcostheta_index[i+1]] + phi_d_2H_Wcostheta_costheta_low[Wcostheta_index[i]:Wcostheta_index[i+1]])/2, \
                    phi_d_2H_Wcostheta_efficiency[0,Wcostheta_index[i]:Wcostheta_index[i+1]], \
                    xerr=(phi_d_2H_Wcostheta_costheta_high[Wcostheta_index[i]:Wcostheta_index[i+1]] - phi_d_2H_Wcostheta_costheta_low[Wcostheta_index[i]:Wcostheta_index[i+1]])/2, \
                    yerr=phi_d_2H_Wcostheta_efficiency_statserr[0,Wcostheta_index[i]:Wcostheta_index[i+1]], \
                    fmt='.', color=W_color_list[i], label=W_label_list[i])
axs[1].set_title(r'$W(\cos \vartheta)$')
axs[1].set_xlabel(r'$\cos\vartheta$')
axs[1].set_ylabel(r'$\epsilon$')
axs[1].set_xlim(-1, 1)
axs[1].set_ylim(0, 0.5)

for i in range(len(Wdecayphi_index)-1):
    axs[2].errorbar((phi_d_2H_Wdecayphi_decayphi_high[Wdecayphi_index[i]:Wdecayphi_index[i+1]] + phi_d_2H_Wdecayphi_decayphi_low[Wdecayphi_index[i]:Wdecayphi_index[i+1]])/2, \
                    phi_d_2H_Wdecayphi_efficiency[0,Wdecayphi_index[i]:Wdecayphi_index[i+1]], \
                    xerr=(phi_d_2H_Wdecayphi_decayphi_high[Wdecayphi_index[i]:Wdecayphi_index[i+1]] - phi_d_2H_Wdecayphi_decayphi_low[Wdecayphi_index[i]:Wdecayphi_index[i+1]])/2, \
                    yerr=phi_d_2H_Wdecayphi_efficiency_statserr[0,Wdecayphi_index[i]:Wdecayphi_index[i+1]], \
                    fmt='.', color=W_color_list[i], label=W_label_list[i])
axs[2].set_title(r'$W(\varphi)$')
axs[2].set_xlabel(r'$\varphi$ (deg)')
axs[2].set_ylabel(r'$\epsilon$')
axs[2].set_xlim(-180, 180)
axs[2].set_ylim(0, 0.5)

for i in range(len(Wpolphi_index)-1):
    axs[3].errorbar((phi_d_2H_Wpolphi_polphi_high[Wpolphi_index[i]:Wpolphi_index[i+1]] + phi_d_2H_Wpolphi_polphi_low[Wpolphi_index[i]:Wpolphi_index[i+1]])/2, \
                    phi_d_2H_Wpolphi_efficiency[0,Wpolphi_index[i]:Wpolphi_index[i+1]], \
                    xerr=(phi_d_2H_Wpolphi_polphi_high[Wpolphi_index[i]:Wpolphi_index[i+1]] - phi_d_2H_Wpolphi_polphi_low[Wpolphi_index[i]:Wpolphi_index[i+1]])/2, \
                    yerr=phi_d_2H_Wpolphi_efficiency_statserr[0,Wpolphi_index[i]:Wpolphi_index[i+1]], \
                    fmt='.', color=W_color_list[i], label=W_label_list[i])
axs[3].set_title(r'$W(\Phi)$')
axs[3].set_xlabel(r'$\Phi$ (deg)')
axs[3].set_ylabel(r'$\epsilon$')
axs[3].set_xlim(-180, 180)
axs[3].set_ylim(0, 0.5)
axs[3].legend(bbox_to_anchor=(1, 0.6), loc='upper right')

for i in range(len(Wpsi_index)-1):
    axs[4].errorbar((phi_d_2H_Wpsi_psi_high[Wpsi_index[i]:Wpsi_index[i+1]] + phi_d_2H_Wpsi_psi_low[Wpsi_index[i]:Wpsi_index[i+1]])/2, \
                    phi_d_2H_Wpsi_efficiency[0,Wpsi_index[i]:Wpsi_index[i+1]], \
                    xerr=(phi_d_2H_Wpsi_psi_high[Wpsi_index[i]:Wpsi_index[i+1]] - phi_d_2H_Wpsi_psi_low[Wpsi_index[i]:Wpsi_index[i+1]])/2, \
                    yerr=phi_d_2H_Wpsi_efficiency_statserr[0,Wpsi_index[i]:Wpsi_index[i+1]], \
                    fmt='.', color=W_color_list[i], label=W_label_list[i])
axs[4].set_title(r'$W(\psi)$')
axs[4].set_xlabel(r'$\psi$ (deg)')
axs[4].set_ylabel(r'$\epsilon$')
axs[4].set_xlim(-180, 180)
axs[4].set_ylim(0, 0.6)

file_pdf.savefig()

################################################################# SYSTEMATICS: DATA YIELD #################################################################################################################

print('Plotting data yield relative difference...')

for k in range(len(cut_name_list)):
    if(cut_name_list[k].find('Simulation parameter') != -1):
        continue

    fig = plt.figure(figsize=(24, 16), dpi=300)
    gs = fig.add_gridspec(2, 6)
    ax1 = fig.add_subplot(gs[0, 0:2])
    ax2 = fig.add_subplot(gs[0, 2:4])
    ax3 = fig.add_subplot(gs[0, 4:6])
    ax4 = fig.add_subplot(gs[1, 1:3])
    ax5 = fig.add_subplot(gs[1, 3:5])
    axs = [ax1, ax2, ax3, ax4, ax5]

    plt.suptitle('Data yield relative difference: '+cut_name_list[k],fontsize=24)
    plt.subplots_adjust(wspace=0.5, hspace=0.5)

    for i in range(1, phi_d_2H_dsdt_yield_data_statserr.shape[0]):
        if (cut_index_list[i] == k):
            for j in range(len(dsdt_index)-1):
                axs[0].scatter( phi_d_2H_dsdt_minust_center[dsdt_index[j]:dsdt_index[j+1]], \
                                (phi_d_2H_dsdt_yield_data[0,dsdt_index[j]:dsdt_index[j+1]]-phi_d_2H_dsdt_yield_data[i,dsdt_index[j]:dsdt_index[j+1]])/(phi_d_2H_dsdt_yield_data[0,dsdt_index[j]:dsdt_index[j+1]]), \
                                label=variation_name_list[i], color=color_list[i])

    for i in range(1, phi_d_2H_Wcostheta_yield_data_statserr.shape[0]):
        if (cut_index_list[i] == k):
            for j in range(len(Wcostheta_index)-1):
                axs[1].scatter( phi_d_2H_Wcostheta_costheta_center[Wcostheta_index[j]:Wcostheta_index[j+1]], \
                                (phi_d_2H_Wcostheta_yield_data[0,Wcostheta_index[j]:Wcostheta_index[j+1]]-phi_d_2H_Wcostheta_yield_data[i,Wcostheta_index[j]:Wcostheta_index[j+1]])/(phi_d_2H_Wcostheta_yield_data[0,Wcostheta_index[j]:Wcostheta_index[j+1]]), \
                                label=variation_name_list[i], color=color_list[i])

    for i in range(1, phi_d_2H_Wdecayphi_yield_data_statserr.shape[0]):
        if (cut_index_list[i] == k):
            for j in range(len(Wdecayphi_index)-1):
                axs[2].scatter( phi_d_2H_Wdecayphi_decayphi_center[Wdecayphi_index[j]:Wdecayphi_index[j+1]], \
                                (phi_d_2H_Wdecayphi_yield_data[0,Wdecayphi_index[j]:Wdecayphi_index[j+1]]-phi_d_2H_Wdecayphi_yield_data[i,Wdecayphi_index[j]:Wdecayphi_index[j+1]])/(phi_d_2H_Wdecayphi_yield_data[0,Wdecayphi_index[j]:Wdecayphi_index[j+1]]), \
                                label=variation_name_list[i], color=color_list[i])

    for i in range(1, phi_d_2H_Wpolphi_yield_data_statserr.shape[0]):
        if (cut_index_list[i] == k):
            for j in range(len(Wpolphi_index)-1):
                axs[3].scatter( phi_d_2H_Wpolphi_polphi_center[Wpolphi_index[j]:Wpolphi_index[j+1]], \
                                (phi_d_2H_Wpolphi_yield_data[0,Wpolphi_index[j]:Wpolphi_index[j+1]]-phi_d_2H_Wpolphi_yield_data[i,Wpolphi_index[j]:Wpolphi_index[j+1]])/(phi_d_2H_Wpolphi_yield_data[0,Wpolphi_index[j]:Wpolphi_index[j+1]]), \
                                label=variation_name_list[i], color=color_list[i])

    for i in range(1, phi_d_2H_Wpsi_yield_data_statserr.shape[0]):
        if (cut_index_list[i] == k):
            for j in range(len(Wpsi_index)-1):
                axs[4].scatter( phi_d_2H_Wpsi_psi_center[Wpsi_index[j]:Wpsi_index[j+1]], \
                                (phi_d_2H_Wpsi_yield_data[0,Wpsi_index[j]:Wpsi_index[j+1]]-phi_d_2H_Wpsi_yield_data[i,Wpsi_index[j]:Wpsi_index[j+1]])/(phi_d_2H_Wpsi_yield_data[0,Wpsi_index[j]:Wpsi_index[j+1]]), \
                                label=variation_name_list[i], color=color_list[i])

    for i, ax in enumerate(axs):
        ax.set_xlim(xmin_list[i], xmax_list[i])
        ax.set_ylim(-0.5, 0.5)
        ax.set_xlabel(xlabel_list[i])
        ax.set_ylabel(r'$(Y^\prime_{\rm data} - Y_{\rm data})/Y_{\rm data}$')
        ax.set_title(title_list[i])
        ax.plot([ax.get_xlim()[0], ax.get_xlim()[1]], [0.1, 0.1], 'r--')
        ax.plot([ax.get_xlim()[0], ax.get_xlim()[1]], [-0.1, -0.1], 'r--')
        if i==2:
            ax.legend()
            legend_without_duplicate_labels(ax)

    file_pdf.savefig()
    plt.close()

################################################################# SYSTEMATICS: SIM YIELD ####################################################################################################################

print('Plotting simulation yield relative difference...')

for k in range(len(cut_name_list)):
    if(cut_name_list[k] == 'Fit max' or cut_name_list[k] == 'Fit width' or cut_name_list[k] == 'Fit background model'):
        continue

    fig = plt.figure(figsize=(24, 16), dpi=300)
    gs = fig.add_gridspec(2, 6)
    ax1 = fig.add_subplot(gs[0, 0:2])
    ax2 = fig.add_subplot(gs[0, 2:4])
    ax3 = fig.add_subplot(gs[0, 4:6])
    ax4 = fig.add_subplot(gs[1, 1:3])
    ax5 = fig.add_subplot(gs[1, 3:5])
    axs = [ax1, ax2, ax3, ax4, ax5]

    plt.suptitle('Simulation yield relative difference: '+cut_name_list[k],fontsize=24)
    plt.subplots_adjust(wspace=0.5, hspace=0.5)

    for i in range(1, phi_d_2H_dsdt_yield_sim_statserr.shape[0]):
        if (cut_index_list[i] == k):
            for j in range(len(dsdt_index)-1):
                axs[0].scatter( phi_d_2H_dsdt_minust_center[dsdt_index[j]:dsdt_index[j+1]], \
                                (phi_d_2H_dsdt_yield_sim[0,dsdt_index[j]:dsdt_index[j+1]]-phi_d_2H_dsdt_yield_sim[i,dsdt_index[j]:dsdt_index[j+1]])/(phi_d_2H_dsdt_yield_sim[0,dsdt_index[j]:dsdt_index[j+1]]), \
                                label=variation_name_list[i], color=color_list[i])

    for i in range(1, phi_d_2H_Wcostheta_yield_sim_statserr.shape[0]):
        if (cut_index_list[i] == k):
            for j in range(len(Wcostheta_index)-1):
                axs[1].scatter( phi_d_2H_Wcostheta_costheta_center[Wcostheta_index[j]:Wcostheta_index[j+1]], \
                                (phi_d_2H_Wcostheta_yield_sim[0,Wcostheta_index[j]:Wcostheta_index[j+1]]-phi_d_2H_Wcostheta_yield_sim[i,Wcostheta_index[j]:Wcostheta_index[j+1]])/(phi_d_2H_Wcostheta_yield_sim[0,Wcostheta_index[j]:Wcostheta_index[j+1]]), \
                                label=variation_name_list[i], color=color_list[i])

    for i in range(1, phi_d_2H_Wdecayphi_yield_sim_statserr.shape[0]):
        if (cut_index_list[i] == k):
            for j in range(len(Wdecayphi_index)-1):
                axs[2].scatter( phi_d_2H_Wdecayphi_decayphi_center[Wdecayphi_index[j]:Wdecayphi_index[j+1]], \
                                (phi_d_2H_Wdecayphi_yield_sim[0,Wdecayphi_index[j]:Wdecayphi_index[j+1]]-phi_d_2H_Wdecayphi_yield_sim[i,Wdecayphi_index[j]:Wdecayphi_index[j+1]])/(phi_d_2H_Wdecayphi_yield_sim[0,Wdecayphi_index[j]:Wdecayphi_index[j+1]]), \
                                label=variation_name_list[i], color=color_list[i])

    for i in range(1, phi_d_2H_Wpolphi_yield_sim_statserr.shape[0]):
        if (cut_index_list[i] == k):
            for j in range(len(Wpolphi_index)-1):
                axs[3].scatter( phi_d_2H_Wpolphi_polphi_center[Wpolphi_index[j]:Wpolphi_index[j+1]], \
                                (phi_d_2H_Wpolphi_yield_sim[0,Wpolphi_index[j]:Wpolphi_index[j+1]]-phi_d_2H_Wpolphi_yield_sim[i,Wpolphi_index[j]:Wpolphi_index[j+1]])/(phi_d_2H_Wpolphi_yield_sim[0,Wpolphi_index[j]:Wpolphi_index[j+1]]), \
                                label=variation_name_list[i], color=color_list[i])

    for i in range(1, phi_d_2H_Wpsi_yield_sim_statserr.shape[0]):
        if (cut_index_list[i] == k):
            for j in range(len(Wpsi_index)-1):
                axs[4].scatter( phi_d_2H_Wpsi_psi_center[Wpsi_index[j]:Wpsi_index[j+1]], \
                                (phi_d_2H_Wpsi_yield_sim[0,Wpsi_index[j]:Wpsi_index[j+1]]-phi_d_2H_Wpsi_yield_sim[i,Wpsi_index[j]:Wpsi_index[j+1]])/(phi_d_2H_Wpsi_yield_sim[0,Wpsi_index[j]:Wpsi_index[j+1]]), \
                                label=variation_name_list[i], color=color_list[i])

    for i, ax in enumerate(axs):
        ax.set_xlim([xmin_list[i], xmax_list[i]])
        ax.set_ylim(-1, 1)
        ax.set_xlabel(xlabel_list[i])
        ax.set_ylabel(r'$(Y^\prime_{\rm sim} - Y_{\rm sim})/Y_{\rm sim}$')
        ax.set_title(title_list[i])
        ax.plot([ax.get_xlim()[0], ax.get_xlim()[1]], [0.1, 0.1],   'r--')
        ax.plot([ax.get_xlim()[0], ax.get_xlim()[1]], [-0.1, -0.1], 'r--')
        if i==2:
            ax.legend()
            legend_without_duplicate_labels(ax)

    file_pdf.savefig()
    plt.close()

################################################################# SYSTEMATICS: THROWN YIELD ####################################################################################################################

print('Plotting thrown yield relative difference...')

for k in range(len(cut_name_list)):
    if(cut_name_list[k].find('Simulation parameter') == -1):
        continue

    fig = plt.figure(figsize=(24, 16), dpi=300)
    gs = fig.add_gridspec(2, 6)
    ax1 = fig.add_subplot(gs[0, 0:2])
    ax2 = fig.add_subplot(gs[0, 2:4])
    ax3 = fig.add_subplot(gs[0, 4:6])
    ax4 = fig.add_subplot(gs[1, 1:3])
    ax5 = fig.add_subplot(gs[1, 3:5])
    axs = [ax1, ax2, ax3, ax4, ax5]

    plt.suptitle('Thrown yield relative difference: '+cut_name_list[k],fontsize=24)
    plt.subplots_adjust(wspace=0.5, hspace=0.5)

    for i in range(1, phi_d_2H_dsdt_yield_tagged_statserr.shape[0]):
        if (cut_index_list[i] == k):
            for j in range(len(dsdt_index)-1):
                axs[0].scatter( phi_d_2H_dsdt_minust_center[dsdt_index[j]:dsdt_index[j+1]], \
                                (phi_d_2H_dsdt_yield_tagged[0,dsdt_index[j]:dsdt_index[j+1]]-phi_d_2H_dsdt_yield_tagged[i,dsdt_index[j]:dsdt_index[j+1]])/(phi_d_2H_dsdt_yield_tagged[0,dsdt_index[j]:dsdt_index[j+1]]), \
                                label=variation_name_list[i], color=color_list[i])

    for i in range(1, phi_d_2H_Wcostheta_yield_tagged_statserr.shape[0]):
        if (cut_index_list[i] == k):
            for j in range(len(Wcostheta_index)-1):
                axs[1].scatter( phi_d_2H_Wcostheta_costheta_center[Wcostheta_index[j]:Wcostheta_index[j+1]], \
                                (phi_d_2H_Wcostheta_yield_tagged[0,Wcostheta_index[j]:Wcostheta_index[j+1]]-phi_d_2H_Wcostheta_yield_tagged[i,Wcostheta_index[j]:Wcostheta_index[j+1]])/(phi_d_2H_Wcostheta_yield_tagged[0,Wcostheta_index[j]:Wcostheta_index[j+1]]), \
                                label=variation_name_list[i], color=color_list[i])

    for i in range(1, phi_d_2H_Wdecayphi_yield_tagged_statserr.shape[0]):
        if (cut_index_list[i] == k):
            for j in range(len(Wdecayphi_index)-1):
                axs[2].scatter( phi_d_2H_Wdecayphi_decayphi_center[Wdecayphi_index[j]:Wdecayphi_index[j+1]], \
                                (phi_d_2H_Wdecayphi_yield_tagged[0,Wdecayphi_index[j]:Wdecayphi_index[j+1]]-phi_d_2H_Wdecayphi_yield_tagged[i,Wdecayphi_index[j]:Wdecayphi_index[j+1]])/(phi_d_2H_Wdecayphi_yield_tagged[0,Wdecayphi_index[j]:Wdecayphi_index[j+1]]), \
                                label=variation_name_list[i], color=color_list[i])

    for i in range(1, phi_d_2H_Wpolphi_yield_tagged_statserr.shape[0]):
        if (cut_index_list[i] == k):
            for j in range(len(Wpolphi_index)-1):
                axs[3].scatter( phi_d_2H_Wpolphi_polphi_center[Wpolphi_index[j]:Wpolphi_index[j+1]], \
                                (phi_d_2H_Wpolphi_yield_tagged[0,Wpolphi_index[j]:Wpolphi_index[j+1]]-phi_d_2H_Wpolphi_yield_tagged[i,Wpolphi_index[j]:Wpolphi_index[j+1]])/(phi_d_2H_Wpolphi_yield_tagged[0,Wpolphi_index[j]:Wpolphi_index[j+1]]), \
                                label=variation_name_list[i], color=color_list[i])

    for i in range(1, phi_d_2H_Wpsi_yield_tagged_statserr.shape[0]):
        if (cut_index_list[i] == k):
            for j in range(len(Wpsi_index)-1):
                axs[4].scatter( phi_d_2H_Wpsi_psi_center[Wpsi_index[j]:Wpsi_index[j+1]], \
                                (phi_d_2H_Wpsi_yield_tagged[0,Wpsi_index[j]:Wpsi_index[j+1]]-phi_d_2H_Wpsi_yield_tagged[i,Wpsi_index[j]:Wpsi_index[j+1]])/(phi_d_2H_Wpsi_yield_tagged[0,Wpsi_index[j]:Wpsi_index[j+1]]), \
                                label=variation_name_list[i], color=color_list[i])

    for i, ax in enumerate(axs):
        ax.set_xlim([xmin_list[i], xmax_list[i]])
        ax.set_ylim(-1, 1)
        ax.set_xlabel(xlabel_list[i])
        ax.set_ylabel(r'$(Y^\prime_{\rm thrown} - Y_{\rm thrown})/Y_{\rm thrown}$')
        ax.set_title(title_list[i])
        ax.plot([ax.get_xlim()[0], ax.get_xlim()[1]], [0.1, 0.1],   'r--')
        ax.plot([ax.get_xlim()[0], ax.get_xlim()[1]], [-0.1, -0.1], 'r--')
        if i==2:
            ax.legend()
            legend_without_duplicate_labels(ax)

    file_pdf.savefig()
    plt.close()

################################################################# SYSTEMATICS: BARLOW SCORE ###############################################################################################################

print('Plotting Barlow score...')

for k in range(len(cut_name_list)):

    fig = plt.figure(figsize=(24, 16), dpi=300)
    gs = fig.add_gridspec(2, 6)
    ax1 = fig.add_subplot(gs[0, 0:2])
    ax2 = fig.add_subplot(gs[0, 2:4])
    ax3 = fig.add_subplot(gs[0, 4:6])
    ax4 = fig.add_subplot(gs[1, 1:3])
    ax5 = fig.add_subplot(gs[1, 3:5])
    axs = [ax1, ax2, ax3, ax4, ax5]

    plt.suptitle('Barlow score: '+cut_name_list[k],fontsize=24)
    plt.subplots_adjust(wspace=0.5, hspace=0.1)

    for i in range(1, phi_d_2H_dsdt_yield_sim_statserr.shape[0]):
        if (cut_index_list[i] == k):
            for j in range(len(dsdt_index)-1):
                axs[0].scatter( phi_d_2H_dsdt_minust_center[dsdt_index[j]:dsdt_index[j+1]], \
                                (phi_d_2H_dsdt_results[0,dsdt_index[j]:dsdt_index[j+1]] - phi_d_2H_dsdt_results[i,dsdt_index[j]:dsdt_index[j+1]])/(1E-40+np.sqrt(np.abs(phi_d_2H_dsdt_results_statserr[0,dsdt_index[j]:dsdt_index[j+1]]**2 - phi_d_2H_dsdt_results_statserr[i,dsdt_index[j]:dsdt_index[j+1]]**2))), \
                                label=variation_name_list[i], color=color_list[i])

    for i in range(1, phi_d_2H_Wcostheta_yield_sim_statserr.shape[0]):
        if (cut_index_list[i] == k):
            for j in range(len(Wcostheta_index)-1):
                axs[1].scatter( phi_d_2H_Wcostheta_costheta_center[Wcostheta_index[j]:Wcostheta_index[j+1]], \
                                (phi_d_2H_Wcostheta_results[0,Wcostheta_index[j]:Wcostheta_index[j+1]] - phi_d_2H_Wcostheta_results[i,Wcostheta_index[j]:Wcostheta_index[j+1]])/(1E-40+np.sqrt(np.abs(phi_d_2H_Wcostheta_results_statserr[0,Wcostheta_index[j]:Wcostheta_index[j+1]]**2 - phi_d_2H_Wcostheta_results_statserr[i,Wcostheta_index[j]:Wcostheta_index[j+1]]**2))), \
                                label=variation_name_list[i], color=color_list[i])

    for i in range(1, phi_d_2H_Wdecayphi_yield_sim_statserr.shape[0]):
        if (cut_index_list[i] == k):
            for j in range(len(Wdecayphi_index)-1):
                axs[2].scatter( phi_d_2H_Wdecayphi_decayphi_center[Wdecayphi_index[j]:Wdecayphi_index[j+1]], \
                                (phi_d_2H_Wdecayphi_results[0,Wdecayphi_index[j]:Wdecayphi_index[j+1]] - phi_d_2H_Wdecayphi_results[i,Wdecayphi_index[j]:Wdecayphi_index[j+1]])/(1E-40+np.sqrt(np.abs(phi_d_2H_Wdecayphi_results_statserr[0,Wdecayphi_index[j]:Wdecayphi_index[j+1]]**2 - phi_d_2H_Wdecayphi_results_statserr[i,Wdecayphi_index[j]:Wdecayphi_index[j+1]]**2))), \
                                label=variation_name_list[i], color=color_list[i])

    for i in range(1, phi_d_2H_Wpolphi_yield_sim_statserr.shape[0]):
        if (cut_index_list[i] == k):
            for j in range(len(Wpolphi_index)-1):
                axs[3].scatter( phi_d_2H_Wpolphi_polphi_center[Wpolphi_index[j]:Wpolphi_index[j+1]], \
                                (phi_d_2H_Wpolphi_results[0,Wpolphi_index[j]:Wpolphi_index[j+1]] - phi_d_2H_Wpolphi_results[i,Wpolphi_index[j]:Wpolphi_index[j+1]])/(1E-40+np.sqrt(np.abs(phi_d_2H_Wpolphi_results_statserr[0,Wpolphi_index[j]:Wpolphi_index[j+1]]**2 - phi_d_2H_Wpolphi_results_statserr[i,Wpolphi_index[j]:Wpolphi_index[j+1]]**2))), \
                                label=variation_name_list[i], color=color_list[i])

    for i in range(1, phi_d_2H_Wpsi_yield_sim_statserr.shape[0]):
        if (cut_index_list[i] == k):
            for j in range(len(Wpsi_index)-1):
                axs[4].scatter( phi_d_2H_Wpsi_psi_center[Wpsi_index[j]:Wpsi_index[j+1]], \
                                (phi_d_2H_Wpsi_results[0,Wpsi_index[j]:Wpsi_index[j+1]]-phi_d_2H_Wpsi_results[i,Wpsi_index[j]:Wpsi_index[j+1]])/(1E-40+np.sqrt(np.abs(phi_d_2H_Wpsi_results_statserr[0,Wpsi_index[j]:Wpsi_index[j+1]]**2 - phi_d_2H_Wpsi_results_statserr[i,Wpsi_index[j]:Wpsi_index[j+1]]**2))), \
                                label=variation_name_list[i], color=color_list[i])

    for i, ax in enumerate(axs):
        ax.set_xlim([xmin_list[i], xmax_list[i]])
        ax.set_ylim(-15, 15)
        ax.set_xlabel(xlabel_list[i])
        ax.set_ylabel(r'Barlow score')
        ax.set_title(title_list[i])
        ax.fill_between([xmin_list[i], xmax_list[i]], [1, 1],   [-1, -1],   color='green',  alpha=0.1)
        ax.fill_between([xmin_list[i], xmax_list[i]], [4, 4],   [1, 1],     color='yellow', alpha=0.1)
        ax.fill_between([xmin_list[i], xmax_list[i]], [15, 15], [4, 4],     color='red',    alpha=0.1)
        ax.fill_between([xmin_list[i], xmax_list[i]], [-1, -1], [-4, -4],   color='yellow', alpha=0.1)
        ax.fill_between([xmin_list[i], xmax_list[i]], [-4, -4], [-15, -15], color='red',    alpha=0.1)
        if i==2:
            ax.legend()
            legend_without_duplicate_labels(ax)

    file_pdf.savefig()
    plt.close()

################################################################# SYSTEMATICS: OBSERVABLE #################################################################################################################

print('Plotting observable relative difference...')

for k in range(len(cut_name_list)):

    fig = plt.figure(figsize=(24, 16), dpi=300)
    gs = fig.add_gridspec(2, 6)
    ax1 = fig.add_subplot(gs[0, 0:2])
    ax2 = fig.add_subplot(gs[0, 2:4])
    ax3 = fig.add_subplot(gs[0, 4:6])
    ax4 = fig.add_subplot(gs[1, 1:3])
    ax5 = fig.add_subplot(gs[1, 3:5])
    axs = [ax1, ax2, ax3, ax4, ax5]

    plt.suptitle('Observable relative difference: '+cut_name_list[k],fontsize=24)
    plt.subplots_adjust(wspace=0.5, hspace=0.5)

    for i in range(1, phi_d_2H_dsdt_yield_sim_statserr.shape[0]):
        if (cut_index_list[i] == k):
            for j in range(len(dsdt_index)-1):
                axs[0].scatter( phi_d_2H_dsdt_minust_center[dsdt_index[j]:dsdt_index[j+1]], \
                                (phi_d_2H_dsdt_results[0,dsdt_index[j]:dsdt_index[j+1]] - phi_d_2H_dsdt_results[i,dsdt_index[j]:dsdt_index[j+1]])/phi_d_2H_dsdt_results[0,dsdt_index[j]:dsdt_index[j+1]], \
                                label=variation_name_list[i], color=color_list[i])

    for i in range(1, phi_d_2H_Wcostheta_yield_sim_statserr.shape[0]):
        if (cut_index_list[i] == k):
            for j in range(len(Wcostheta_index)-1):
                axs[1].scatter( phi_d_2H_Wcostheta_costheta_center[Wcostheta_index[j]:Wcostheta_index[j+1]], \
                                (phi_d_2H_Wcostheta_results[0,Wcostheta_index[j]:Wcostheta_index[j+1]] - phi_d_2H_Wcostheta_results[i,Wcostheta_index[j]:Wcostheta_index[j+1]])/phi_d_2H_Wcostheta_results[0,Wcostheta_index[j]:Wcostheta_index[j+1]], \
                                label=variation_name_list[i], color=color_list[i])

    for i in range(1, phi_d_2H_Wdecayphi_yield_sim_statserr.shape[0]):
        if (cut_index_list[i] == k):
            for j in range(len(Wdecayphi_index)-1):
                axs[2].scatter( phi_d_2H_Wdecayphi_decayphi_center[Wdecayphi_index[j]:Wdecayphi_index[j+1]], \
                                (phi_d_2H_Wdecayphi_results[0,Wdecayphi_index[j]:Wdecayphi_index[j+1]] - phi_d_2H_Wdecayphi_results[i,Wdecayphi_index[j]:Wdecayphi_index[j+1]])/phi_d_2H_Wdecayphi_results[0,Wdecayphi_index[j]:Wdecayphi_index[j+1]], \
                                label=variation_name_list[i], color=color_list[i])

    for i in range(1, phi_d_2H_Wpolphi_yield_sim_statserr.shape[0]):
        if (cut_index_list[i] == k):
            for j in range(len(Wpolphi_index)-1):
                axs[3].scatter( phi_d_2H_Wpolphi_polphi_center[Wpolphi_index[j]:Wpolphi_index[j+1]], \
                                (phi_d_2H_Wpolphi_results[0,Wpolphi_index[j]:Wpolphi_index[j+1]] - phi_d_2H_Wpolphi_results[i,Wpolphi_index[j]:Wpolphi_index[j+1]])/phi_d_2H_Wpolphi_results[0,Wpolphi_index[j]:Wpolphi_index[j+1]], \
                                label=variation_name_list[i], color=color_list[i])

    for i in range(1, phi_d_2H_Wpsi_yield_sim_statserr.shape[0]):
        if (cut_index_list[i] == k):
            for j in range(len(Wpsi_index)-1):
                axs[4].scatter( phi_d_2H_Wpsi_psi_center[Wpsi_index[j]:Wpsi_index[j+1]], \
                                (phi_d_2H_Wpsi_results[0,Wpsi_index[j]:Wpsi_index[j+1]]-phi_d_2H_Wpsi_results[i,Wpsi_index[j]:Wpsi_index[j+1]])/phi_d_2H_Wpsi_results[0,Wpsi_index[j]:Wpsi_index[j+1]], \
                                label=variation_name_list[i], color=color_list[i])

    for i, ax in enumerate(axs):
        ax.set_xlim([xmin_list[i], xmax_list[i]])
        ax.set_ylim(-0.2, 0.2)
        ax.set_xlabel(xlabel_list[i])
        ax.set_ylabel(r'$(\sigma^\prime - \sigma)/\sigma$')
        ax.set_title(title_list[i])
        if i==2:
            ax.legend()
            legend_without_duplicate_labels(ax)

    file_pdf.savefig()
    plt.close()

################################################################# SYSTEMATICS: DSDT #############################################################################################################

print('Plotting systematic uncertainties...')

# Point-to-point
phi_d_2H_dsdt_results_p2perr    = np.zeros((len(p2p_list), np.shape(phi_d_2H_dsdt_results)[1]))
for i in range(phi_d_2H_dsdt_yield_sim_statserr.shape[0]):
    if (i == 0):
        temp_array = np.zeros_like(phi_d_2H_dsdt_results[0])
        continue
    temp_array = np.vstack((temp_array, (phi_d_2H_dsdt_results[0] - phi_d_2H_dsdt_results[i])/phi_d_2H_dsdt_results[0]))
    if (i == phi_d_2H_dsdt_yield_sim_statserr.shape[0]-1 or cut_index_list[i+1] != cut_index_list[i]):
        temp_std = np.std(temp_array, axis=0)
        phi_d_2H_dsdt_results_p2perr[cut_class_list[cut_index_list[i]]] += temp_std**2
        temp_array = np.zeros_like(phi_d_2H_dsdt_results[0])

for i in range(len(p2p_list)):
    if p2p_list[i] == 'Simulation model':  # No systematics assigned due to the simulation model
        continue
    phi_d_2H_dsdt_results_p2perr[0] += phi_d_2H_dsdt_results_p2perr[i]
phi_d_2H_dsdt_results_p2perr        = phi_d_2H_dsdt_results[0]*np.sqrt(phi_d_2H_dsdt_results_p2perr)

# Normalization
phi_d_2H_dsdt_results_normerr       = 0
phi_d_2H_dsdt_results_normerr       += 0.0656**2            # Tracking efficiency
phi_d_2H_dsdt_results_normerr       += 0.05**2              # Flux measurement
phi_d_2H_dsdt_results_normerr       += 0.002**2 + 0.005**2  # Target properties
phi_d_2H_dsdt_results_normerr       += 0.01**2              # Branching ratio
phi_d_2H_dsdt_results_normerr       = np.sqrt(phi_d_2H_dsdt_results_normerr)
phi_d_2H_dsdt_results_normerr       = phi_d_2H_dsdt_results[0]*phi_d_2H_dsdt_results_normerr

# Total
phi_d_2H_dsdt_results_systerr      = np.sqrt(phi_d_2H_dsdt_results_p2perr[0]**2 + phi_d_2H_dsdt_results_normerr**2)
phi_d_2H_dsdt_results_totalerr     = np.sqrt(phi_d_2H_dsdt_results_statserr[0]**2 + phi_d_2H_dsdt_results_systerr**2)

# Plot
fig = plt.figure(figsize=(24, 8), dpi=300)
gs = fig.add_gridspec(1, 3, wspace=0)
axs = gs.subplots(sharey=True)

for j in range(len(dsdt_index)-1):
    axs[j].scatter(phi_d_2H_dsdt_minust_center[dsdt_index[j]:dsdt_index[j+1]], phi_d_2H_dsdt_results_systerr[dsdt_index[j]:dsdt_index[j+1]]/phi_d_2H_dsdt_results[0, dsdt_index[j]:dsdt_index[j+1]], label='Total systematic', color='black', s=100)
    axs[j].text(0.05, 0.95, dsdt_label_list[j], transform=axs[j].transAxes, fontsize=16, verticalalignment='top')
    axs[j].set_xlim(0, 2)
    axs[j].set_xlabel(xlabel_list[0], fontsize=16)
    axs[j].tick_params(labelsize=16)
axs[0].set_ylabel(r'$\delta\sigma/\sigma$', fontsize=16)
axs[0].set_ylim(0, 0.2)
axs[0].set_xticks(np.arange(0, 1.8, 0.25))
axs[1].set_xticks(np.arange(0, 1.8, 0.25))
axs[2].set_xticks(np.arange(0, 2.1, 0.25))

for i in range(len(p2p_list)):
    if p2p_list[i] == 'Simulation model':
        continue
    for j in range(len(dsdt_index)-1):
        axs[j].scatter(phi_d_2H_dsdt_minust_center[dsdt_index[j]:dsdt_index[j+1]], phi_d_2H_dsdt_results_p2perr[i, dsdt_index[j]:dsdt_index[j+1]]/phi_d_2H_dsdt_results[0, dsdt_index[j]:dsdt_index[j+1]], label=p2p_list[i], marker='x' if i!=0 else 's', s=100)

plt.suptitle('Systematic uncertainties: ' + title_list[0], fontsize=20)
plt.legend()
file_pdf.savefig()
plt.close()

################################################################# SYSTEMATICS: W(COSTHETA) #############################################################################################################

# Point-to-point
phi_d_2H_Wcostheta_results_p2perr    = np.zeros((len(p2p_list), np.shape(phi_d_2H_Wcostheta_results)[1]))
for i in range(phi_d_2H_Wcostheta_yield_sim_statserr.shape[0]):
    if (i == 0):
        temp_array = np.zeros_like(phi_d_2H_Wcostheta_results[0])
        continue
    temp_array = np.vstack((temp_array, (phi_d_2H_Wcostheta_results[0] - phi_d_2H_Wcostheta_results[i])/phi_d_2H_Wcostheta_results[0]))
    if (i == phi_d_2H_Wcostheta_yield_sim_statserr.shape[0]-1 or cut_index_list[i+1] != cut_index_list[i]):
        temp_std = np.std(temp_array, axis=0)
        phi_d_2H_Wcostheta_results_p2perr[cut_class_list[cut_index_list[i]]] += temp_std**2
        temp_array = np.zeros_like(phi_d_2H_Wcostheta_results[0])

phi_d_2H_Wcostheta_results_p2perr[0]    = np.sum(phi_d_2H_Wcostheta_results_p2perr[1:], axis=0)
phi_d_2H_Wcostheta_results_p2perr       = phi_d_2H_Wcostheta_results[0]*np.sqrt(phi_d_2H_Wcostheta_results_p2perr)

# Normalization
phi_d_2H_Wcostheta_results_normerr       = 0

# Total
phi_d_2H_Wcostheta_results_systerr       = np.sqrt(phi_d_2H_Wcostheta_results_p2perr[0]**2 + phi_d_2H_Wcostheta_results_normerr**2)
phi_d_2H_Wcostheta_results_totalerr      = np.sqrt(phi_d_2H_Wcostheta_results_statserr[0]**2 + phi_d_2H_Wcostheta_results_systerr**2)

# Plot
fig = plt.figure(figsize=(24, 16), dpi=300)
gs = fig.add_gridspec(2, 3, wspace=0, hspace=0)
axs = gs.subplots(sharey=True).flatten()

for j in range(len(Wcostheta_index)-1):
    axs[j].scatter(phi_d_2H_Wcostheta_costheta_center[Wcostheta_index[j]:Wcostheta_index[j+1]], phi_d_2H_Wcostheta_results_systerr[Wcostheta_index[j]:Wcostheta_index[j+1]]/phi_d_2H_Wcostheta_results[0, Wcostheta_index[j]:Wcostheta_index[j+1]], s=100, label='Total systematic', color='black')
    axs[j].text(0.05, 0.95, W_label_list[j], transform=axs[j].transAxes, fontsize=16, verticalalignment='top')
    axs[j].set_xlim(-1, 1)
    axs[j].tick_params(labelsize=16)

for i in range(len(p2p_list)):
    if i==0:
        continue
    for j in range(len(Wcostheta_index)-1):
        axs[j].scatter(phi_d_2H_Wcostheta_costheta_center[Wcostheta_index[j]:Wcostheta_index[j+1]], phi_d_2H_Wcostheta_results_p2perr[i, Wcostheta_index[j]:Wcostheta_index[j+1]]/phi_d_2H_Wcostheta_results[0, Wcostheta_index[j]:Wcostheta_index[j+1]], s=100, label=p2p_list[i], marker='x' if i!=0 else 's')

axs[0].set_ylabel(r'$\delta W/W$', fontsize=16)
axs[3].set_ylabel(r'$\delta W/W$', fontsize=16)
axs[0].set_ylim(0, 0.2)
axs[0].set_xticks([])
axs[1].set_xticks([])
axs[2].set_xticks([])
axs[3].set_xticks(np.arange(-1, 1, 0.25))
axs[4].set_xticks(np.arange(-1, 1, 0.25))
axs[5].set_xticks(np.arange(-1, 1.1, 0.25))
axs[3].set_xlabel(xlabel_list[1], fontsize=16)
axs[4].set_xlabel(xlabel_list[1], fontsize=16)
axs[5].set_xlabel(xlabel_list[1], fontsize=16)
axs[2].legend()

plt.suptitle('Systematic uncertainties: ' + title_list[1], fontsize=20)
file_pdf.savefig()
plt.close()

################################################################# SYSTEMATICS: W(DECAYPHI) #############################################################################################################

# Point-to-point
phi_d_2H_Wdecayphi_results_p2perr    = np.zeros((len(p2p_list), np.shape(phi_d_2H_Wdecayphi_results)[1]))
for i in range(phi_d_2H_Wdecayphi_yield_sim_statserr.shape[0]):
    if (i == 0):
        temp_array = np.zeros_like(phi_d_2H_Wdecayphi_results[0])
        continue
    temp_array = np.vstack((temp_array, (phi_d_2H_Wdecayphi_results[0] - phi_d_2H_Wdecayphi_results[i])/phi_d_2H_Wdecayphi_results[0]))
    if (i == phi_d_2H_Wdecayphi_yield_sim_statserr.shape[0]-1 or cut_index_list[i+1] != cut_index_list[i]):
        temp_std = np.std(temp_array, axis=0)
        phi_d_2H_Wdecayphi_results_p2perr[cut_class_list[cut_index_list[i]]] += temp_std**2
        temp_array = np.zeros_like(phi_d_2H_Wdecayphi_results[0])

phi_d_2H_Wdecayphi_results_p2perr[0]    = np.sum(phi_d_2H_Wdecayphi_results_p2perr[1:], axis=0)
phi_d_2H_Wdecayphi_results_p2perr       = phi_d_2H_Wdecayphi_results[0]*np.sqrt(phi_d_2H_Wdecayphi_results_p2perr)

# Normalization
phi_d_2H_Wdecayphi_results_normerr      = 0

# Total
phi_d_2H_Wdecayphi_results_systerr      = np.sqrt(phi_d_2H_Wdecayphi_results_p2perr[0]**2 + phi_d_2H_Wdecayphi_results_normerr**2)
phi_d_2H_Wdecayphi_results_totalerr     = np.sqrt(phi_d_2H_Wdecayphi_results_statserr[0]**2 + phi_d_2H_Wdecayphi_results_systerr**2)

# Plot
fig = plt.figure(figsize=(24, 16), dpi=300)
gs = fig.add_gridspec(2, 3, wspace=0, hspace=0)
axs = gs.subplots(sharey=True).flatten()

for j in range(len(Wdecayphi_index)-1):
    axs[j].scatter(phi_d_2H_Wdecayphi_decayphi_center[Wdecayphi_index[j]:Wdecayphi_index[j+1]], phi_d_2H_Wdecayphi_results_systerr[Wdecayphi_index[j]:Wdecayphi_index[j+1]]/phi_d_2H_Wdecayphi_results[0, Wdecayphi_index[j]:Wdecayphi_index[j+1]], label='Total systematic', color='black', s=100)
    axs[j].text(0.05, 0.95, W_label_list[j], transform=axs[j].transAxes, fontsize=16, verticalalignment='top')
    axs[j].set_xlim(-180, 180)
    axs[j].tick_params(labelsize=16)

for i in range(len(p2p_list)):
    if i==0:
        continue
    for j in range(len(Wdecayphi_index)-1):
        axs[j].scatter(phi_d_2H_Wdecayphi_decayphi_center[Wdecayphi_index[j]:Wdecayphi_index[j+1]], phi_d_2H_Wdecayphi_results_p2perr[i, Wdecayphi_index[j]:Wdecayphi_index[j+1]]/phi_d_2H_Wdecayphi_results[0, Wdecayphi_index[j]:Wdecayphi_index[j+1]], label=p2p_list[i], marker='x' if i!=0 else 's', s=100)

axs[0].set_ylabel(r'$\delta W/W$', fontsize=16)
axs[3].set_ylabel(r'$\delta W/W$', fontsize=16)
axs[0].set_ylim(0, 0.2)
axs[0].set_xticks([])
axs[1].set_xticks([])
axs[2].set_xticks([])
axs[3].set_xticks(np.arange(-180, 180, 45))
axs[4].set_xticks(np.arange(-180, 180, 45))
axs[5].set_xticks(np.arange(-180, 181, 45))
axs[3].set_xlabel(xlabel_list[2], fontsize=16)
axs[4].set_xlabel(xlabel_list[2], fontsize=16)
axs[5].set_xlabel(xlabel_list[2], fontsize=16)
axs[2].legend()

plt.suptitle('Systematic uncertainties: ' + title_list[2], fontsize=20)
file_pdf.savefig()
plt.close()

################################################################# SYSTEMATICS: W(POLPHI) #############################################################################################################

# Point-to-point
phi_d_2H_Wpolphi_results_p2perr    = np.zeros((len(p2p_list), np.shape(phi_d_2H_Wpolphi_results)[1]))
for i in range(phi_d_2H_Wpolphi_yield_sim_statserr.shape[0]):
    if (i == 0):
        temp_array = np.zeros_like(phi_d_2H_Wpolphi_results[0])
        continue
    temp_array = np.vstack((temp_array, (phi_d_2H_Wpolphi_results[0] - phi_d_2H_Wpolphi_results[i])/phi_d_2H_Wpolphi_results[0]))
    if (i == phi_d_2H_Wpolphi_yield_sim_statserr.shape[0]-1 or cut_index_list[i+1] != cut_index_list[i]):
        temp_std = np.std(temp_array, axis=0)
        phi_d_2H_Wpolphi_results_p2perr[cut_class_list[cut_index_list[i]]] += temp_std**2
        temp_array = np.zeros_like(phi_d_2H_Wpolphi_results[0])

phi_d_2H_Wpolphi_results_p2perr[0]    = np.sum(phi_d_2H_Wpolphi_results_p2perr[1:], axis=0)
phi_d_2H_Wpolphi_results_p2perr       = phi_d_2H_Wpolphi_results[0]*np.sqrt(phi_d_2H_Wpolphi_results_p2perr)

# Normalization
phi_d_2H_Wpolphi_results_normerr      = 0

# Total
phi_d_2H_Wpolphi_results_systerr     = np.sqrt(phi_d_2H_Wpolphi_results_p2perr[0]**2 + phi_d_2H_Wpolphi_results_normerr**2)
phi_d_2H_Wpolphi_results_totalerr    = np.sqrt(phi_d_2H_Wpolphi_results_statserr[0]**2 + phi_d_2H_Wpolphi_results_systerr**2)

# Plot
fig = plt.figure(figsize=(24, 16), dpi=300)
gs = fig.add_gridspec(2, 3, wspace=0, hspace=0)
axs = gs.subplots(sharey=True).flatten()

for j in range(len(Wpolphi_index)-1):
    axs[j].scatter(phi_d_2H_Wpolphi_polphi_center[Wpolphi_index[j]:Wpolphi_index[j+1]], phi_d_2H_Wpolphi_results_systerr[Wpolphi_index[j]:Wpolphi_index[j+1]]/phi_d_2H_Wpolphi_results[0, Wpolphi_index[j]:Wpolphi_index[j+1]], label='Total systematic', color='black', s=100)
    axs[j].text(0.05, 0.95, W_label_list[j], transform=axs[j].transAxes, fontsize=16, verticalalignment='top')
    axs[j].set_xlim(-180, 180)
    axs[j].tick_params(labelsize=16)

for i in range(len(p2p_list)):
    if i==0:
        continue
    for j in range(len(Wpolphi_index)-1):
        axs[j].scatter(phi_d_2H_Wpolphi_polphi_center[Wpolphi_index[j]:Wpolphi_index[j+1]], phi_d_2H_Wpolphi_results_p2perr[i, Wpolphi_index[j]:Wpolphi_index[j+1]]/phi_d_2H_Wpolphi_results[0, Wpolphi_index[j]:Wpolphi_index[j+1]], label=p2p_list[i], marker='x' if i!=0 else 's', s=100)

axs[0].set_ylabel(r'$\delta W/W$', fontsize=16)
axs[3].set_ylabel(r'$\delta W/W$', fontsize=16)
axs[0].set_ylim(0, 0.2)
axs[0].set_xticks([])
axs[1].set_xticks([])
axs[2].set_xticks([])
axs[3].set_xticks(np.arange(-180, 180, 45))
axs[4].set_xticks(np.arange(-180, 180, 45))
axs[5].set_xticks(np.arange(-180, 181, 45))
axs[3].set_xlabel(xlabel_list[3], fontsize=16)
axs[4].set_xlabel(xlabel_list[3], fontsize=16)
axs[5].set_xlabel(xlabel_list[3], fontsize=16)
axs[2].legend()

plt.suptitle('Systematic uncertainties: ' + title_list[3], fontsize=20)
file_pdf.savefig()
plt.close()

################################################################# SYSTEMATICS: W(PSI) #############################################################################################################

# Point-to-point
phi_d_2H_Wpsi_results_p2perr    = np.zeros((len(p2p_list), np.shape(phi_d_2H_Wpsi_results)[1]))
for i in range(phi_d_2H_Wpsi_yield_sim_statserr.shape[0]):
    if (i == 0):
        temp_array = np.zeros_like(phi_d_2H_Wpsi_results[0])
        continue
    temp_array = np.vstack((temp_array, (phi_d_2H_Wpsi_results[0] - phi_d_2H_Wpsi_results[i])/phi_d_2H_Wpsi_results[0]))
    if (i == phi_d_2H_Wpsi_yield_sim_statserr.shape[0]-1 or cut_index_list[i+1] != cut_index_list[i]):
        temp_std = np.std(temp_array, axis=0)
        phi_d_2H_Wpsi_results_p2perr[cut_class_list[cut_index_list[i]]] += temp_std**2
        temp_array = np.zeros_like(phi_d_2H_Wpsi_results[0])

phi_d_2H_Wpsi_results_p2perr[0]    = np.sum(phi_d_2H_Wpsi_results_p2perr[1:], axis=0)
phi_d_2H_Wpsi_results_p2perr       = phi_d_2H_Wpsi_results[0]*np.sqrt(phi_d_2H_Wpsi_results_p2perr)

# Normalization
phi_d_2H_Wpsi_results_normerr      = 0

# Total
phi_d_2H_Wpsi_results_systerr      = np.sqrt(phi_d_2H_Wpsi_results_p2perr[0]**2 + phi_d_2H_Wpsi_results_normerr**2)
phi_d_2H_Wpsi_results_totalerr     = np.sqrt(phi_d_2H_Wpsi_results_statserr[0]**2 + phi_d_2H_Wpsi_results_systerr**2)

# Plot
fig = plt.figure(figsize=(24, 16), dpi=300)
gs = fig.add_gridspec(2, 3, wspace=0, hspace=0)
axs = gs.subplots(sharey=True).flatten()

for j in range(len(Wpsi_index)-1):
    axs[j].scatter(phi_d_2H_Wpsi_psi_center[Wpsi_index[j]:Wpsi_index[j+1]], phi_d_2H_Wpsi_results_systerr[Wpsi_index[j]:Wpsi_index[j+1]]/phi_d_2H_Wpsi_results[0, Wpsi_index[j]:Wpsi_index[j+1]], label='Total systematic', color='black', s=100)
    axs[j].text(0.05, 0.95, W_label_list[j], transform=axs[j].transAxes, fontsize=16, verticalalignment='top')
    axs[j].set_xlim(-180, 180)
    axs[j].tick_params(labelsize=16)

for i in range(len(p2p_list)):
    if i==0:
        continue
    for j in range(len(Wpsi_index)-1):
        axs[j].scatter(phi_d_2H_Wpsi_psi_center[Wpsi_index[j]:Wpsi_index[j+1]], phi_d_2H_Wpsi_results_p2perr[i, Wpsi_index[j]:Wpsi_index[j+1]]/phi_d_2H_Wpsi_results[0, Wpsi_index[j]:Wpsi_index[j+1]], label=p2p_list[i], marker='x' if i!=0 else 's', s=100)

axs[0].set_ylabel(r'$\delta W/W$', fontsize=16)
axs[3].set_ylabel(r'$\delta W/W$', fontsize=16)
axs[0].set_ylim(0, 0.2)
axs[0].set_xticks([])
axs[1].set_xticks([])
axs[2].set_xticks([])
axs[3].set_xticks(np.arange(-180, 180, 45))
axs[4].set_xticks(np.arange(-180, 180, 45))
axs[5].set_xticks(np.arange(-180, 181, 45))
axs[3].set_xlabel(xlabel_list[4], fontsize=16)
axs[4].set_xlabel(xlabel_list[4], fontsize=16)
axs[5].set_xlabel(xlabel_list[4], fontsize=16)
axs[2].legend()

plt.suptitle('Systematic uncertainties: ' + title_list[4], fontsize=20)
file_pdf.savefig()
plt.close()

################################################################# OBSERVABLE RESULTS: DIFFERENTIAL CROSS SECTIONS #############################################################################################################

print('Plotting dsdt results...')

# Data points from CLAS
phi_d_2H_dsdt_clas_minust_low           = np.array([0.350, 0.375, 0.400, 0.425, 0.450, 0.500, 0.550, 0.600, 0.700, 0.800, 1.000, 1.200, 1.400])
phi_d_2H_dsdt_clas_minust_high          = np.array([0.375, 0.400, 0.425, 0.450, 0.500, 0.550, 0.600, 0.700, 0.800, 1.000, 1.200, 1.400, 2.000])
phi_d_2H_dsdt_clas_minust_middle        = (phi_d_2H_dsdt_clas_minust_high + phi_d_2H_dsdt_clas_minust_low) / 2
phi_d_2H_dsdt_clas_minust_size          = (phi_d_2H_dsdt_clas_minust_high - phi_d_2H_dsdt_clas_minust_low) / 2
phi_d_2H_dsdt_clas_minust_center        = np.array([0.360, 0.385, 0.410, 0.435, 0.474, 0.524, 0.574, 0.646, 0.746, 0.888, 1.091, 1.292, 1.637])
phi_d_2H_dsdt_clas_results_16           = np.array([10.21, 8.85, 7.32, 6.16, 4.73, 3.52, 2.66, 2.17, 1.40, 0.94, 0.57, 0.28, 0.19])
phi_d_2H_dsdt_clas_results_16_statserr  = np.array([0.82, 0.75, 0.59, 0.55, 0.34, 0.28, 0.24, 0.15, 0.12, 0.07, 0.06, 0.05, 0.02])
phi_d_2H_dsdt_clas_results_16_systerr   = np.array([1.70, 1.11, 0.94, 0.81, 0.60, 0.51, 0.38, 0.26, 0.16, 0.11, 0.07, 0.04, 0.03])
phi_d_2H_dsdt_clas_results_16_totalerr  = np.sqrt(phi_d_2H_dsdt_clas_results_16_statserr**2 + phi_d_2H_dsdt_clas_results_16_systerr**2)
phi_d_2H_dsdt_clas_results_26           = np.array([8.63, 6.80, 4.57, 5.76, 3.99, 3.59, 2.11, 1.83, 1.32, 0.96, 0.57, 0.36, 0.15])
phi_d_2H_dsdt_clas_results_26_statserr  = np.array([0.80, 0.69, 0.53, 0.56, 0.33, 0.29, 0.22, 0.14, 0.12, 0.07, 0.05, 0.04, 0.02])
phi_d_2H_dsdt_clas_results_26_systerr   = np.array([1.04, 1.07, 0.74, 0.65, 0.55, 0.55, 0.28, 0.24, 0.20, 0.11, 0.06, 0.05, 0.02])
phi_d_2H_dsdt_clas_results_26_totalerr  = np.sqrt(phi_d_2H_dsdt_clas_results_26_statserr**2 + phi_d_2H_dsdt_clas_results_26_systerr**2)

# Data points from LEPS
phi_d_2H_dsdt_leps_minust_low           = np.array([0.38, 0.36, 0.34, 0.32, 0.30, 0.28, 0.26, 0.24, 0.22, 0.20, 0.18, 0.16, 0.14, 0.12, 0.10])
phi_d_2H_dsdt_leps_minust_high          = np.array([0.4, 0.38, 0.36, 0.34, 0.32, 0.30, 0.28, 0.26, 0.24, 0.22, 0.20, 0.18, 0.16, 0.14, 0.12])
phi_d_2H_dsdt_leps_minust_center        = (phi_d_2H_dsdt_leps_minust_high + phi_d_2H_dsdt_leps_minust_low) / 2
phi_d_2H_dsdt_leps_minust_width         = (phi_d_2H_dsdt_leps_minust_high - phi_d_2H_dsdt_leps_minust_low) / 2
phi_d_2H_dsdt_leps_results_157          = np.array([0.0005, 0.004, 0.0087, 0.0068, 0.0238, 0.0317, 0.0567, 0.0722, 0.092, 0.1186, 0.1749, 0.2033, 0.2544, 0.3101, 0.3396])*1000
phi_d_2H_dsdt_leps_results_157_statserr = np.array([0.0005, 0.002, 0.0035, 0.0028, 0.007, 0.0076, 0.0102, 0.0118, 0.0142, 0.0137, 0.0159, 0.0148, 0.0166, 0.0152, 0.0143])*1000

# Theoretical calculation examples
phi_d_2H_dsdt_theory_minust             = np.loadtxt('/work/halld2/home/boyu/src_analysis/plot/vm_d/theory/output/nominal/E_8.3_sigma_30.0_b_10.0.txt')[:,0]
phi_d_2H_dsdt_theory_10mb               = np.loadtxt('/work/halld2/home/boyu/src_analysis/plot/vm_d/theory/output/nominal/E_8.3_sigma_20.0_b_6.0.txt')[:,2]
phi_d_2H_dsdt_theory_30mb               = np.loadtxt('/work/halld2/home/boyu/src_analysis/plot/vm_d/theory/output/nominal/E_8.3_sigma_30.0_b_8.0.txt')[:,2]

# ds/dt world data
fig = plt.figure(figsize=(8, 6), dpi=300)
plt.errorbar(phi_d_2H_dsdt_minust_center[dsdt_index[0]:dsdt_index[1]],  phi_d_2H_dsdt_results[0, dsdt_index[0]:dsdt_index[1]],  xerr=phi_d_2H_dsdt_minust_width[dsdt_index[0]:dsdt_index[1]],   yerr=phi_d_2H_dsdt_results_statserr[0, dsdt_index[0]:dsdt_index[1]],    fmt='b.', capsize=1.5, capthick=1, label='This work (5.8-7.8 GeV)')
plt.errorbar(phi_d_2H_dsdt_minust_center[dsdt_index[1]:dsdt_index[2]],  phi_d_2H_dsdt_results[0, dsdt_index[1]:dsdt_index[2]],  xerr=phi_d_2H_dsdt_minust_width[dsdt_index[1]:dsdt_index[2]],   yerr=phi_d_2H_dsdt_results_statserr[0, dsdt_index[1]:dsdt_index[2]],    fmt='k.', capsize=1.5, capthick=1, label='This work (7.8-8.8 GeV)')
plt.errorbar(phi_d_2H_dsdt_minust_center[dsdt_index[2]:dsdt_index[3]],  phi_d_2H_dsdt_results[0, dsdt_index[2]:dsdt_index[3]],  xerr=phi_d_2H_dsdt_minust_width[dsdt_index[2]:dsdt_index[3]],   yerr=phi_d_2H_dsdt_results_statserr[0, dsdt_index[2]:dsdt_index[3]],    fmt='r.', capsize=1.5, capthick=1, label='This work (8.8-10.8 GeV)')
plt.errorbar(phi_d_2H_dsdt_minust_center[dsdt_index[1]:dsdt_index[2]],  phi_d_2H_dsdt_results[0, dsdt_index[1]:dsdt_index[2]],  xerr=phi_d_2H_dsdt_minust_width[dsdt_index[1]:dsdt_index[2]],   yerr=phi_d_2H_dsdt_results_totalerr[dsdt_index[1]:dsdt_index[2]],       fmt='k.')
plt.errorbar(phi_d_2H_dsdt_minust_center[dsdt_index[2]:dsdt_index[3]],  phi_d_2H_dsdt_results[0, dsdt_index[2]:dsdt_index[3]],  xerr=phi_d_2H_dsdt_minust_width[dsdt_index[2]:dsdt_index[3]],   yerr=phi_d_2H_dsdt_results_totalerr[dsdt_index[2]:dsdt_index[3]],       fmt='r.')

plt.errorbar(phi_d_2H_dsdt_clas_minust_center,                          phi_d_2H_dsdt_clas_results_16,                                                                                          yerr=phi_d_2H_dsdt_clas_results_16_statserr,                            fmt='ys', markerfacecolor='white', markersize=3, capsize=2, capthick=1)
plt.errorbar(phi_d_2H_dsdt_clas_minust_center,                          phi_d_2H_dsdt_clas_results_16,                                                                                          yerr=phi_d_2H_dsdt_clas_results_16_totalerr,                            fmt='ys', markerfacecolor='white', markersize=3, label='CLAS (1.6-2.6 GeV)')
plt.errorbar(phi_d_2H_dsdt_clas_minust_center,                          phi_d_2H_dsdt_clas_results_26,                                                                                          yerr=phi_d_2H_dsdt_clas_results_26_statserr,                            fmt='cs', markerfacecolor='white', markersize=3, capsize=2, capthick=1)
plt.errorbar(phi_d_2H_dsdt_clas_minust_center,                          phi_d_2H_dsdt_clas_results_26,                                                                                          yerr=phi_d_2H_dsdt_clas_results_26_totalerr,                            fmt='cs', markerfacecolor='white', markersize=3, label='CLAS (2.6-3.6 GeV)')
plt.errorbar(phi_d_2H_dsdt_leps_minust_center,                          phi_d_2H_dsdt_leps_results_157,                                                                                         yerr=phi_d_2H_dsdt_leps_results_157_statserr,                           fmt='g^', markerfacecolor='white', markersize=3, capsize=2, capthick=1)
plt.errorbar(phi_d_2H_dsdt_leps_minust_center,                          phi_d_2H_dsdt_leps_results_157,                                                                                         yerr=phi_d_2H_dsdt_leps_results_157_statserr,                           fmt='g^', markerfacecolor='white', markersize=3, label='LEPS (1.57-2.37 GeV)')

plt.xlabel(r'$-t\ [GeV^2/c]$', fontsize=16)
plt.ylabel(r'$d\sigma/dt\ [nb/(GeV^2/c)]$', fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.xlim(0, 2)
plt.ylim(1e-1, 1e3)
plt.yscale('log')
plt.legend()
file_pdf.savefig()
plt.close()

# ds/dt compare with theory
fig = plt.figure(figsize=(8, 6), dpi=300)
plt.errorbar(phi_d_2H_dsdt_minust_center[dsdt_index[0]:dsdt_index[1]],  phi_d_2H_dsdt_results[0, dsdt_index[0]:dsdt_index[1]],  xerr=phi_d_2H_dsdt_minust_width[dsdt_index[0]:dsdt_index[1]],   yerr=phi_d_2H_dsdt_results_statserr[0, dsdt_index[0]:dsdt_index[1]],    fmt='b.', capsize=1.5, capthick=1, label='This work (5.8-7.8 GeV)')
plt.errorbar(phi_d_2H_dsdt_minust_center[dsdt_index[1]:dsdt_index[2]],  phi_d_2H_dsdt_results[0, dsdt_index[1]:dsdt_index[2]],  xerr=phi_d_2H_dsdt_minust_width[dsdt_index[1]:dsdt_index[2]],   yerr=phi_d_2H_dsdt_results_statserr[0, dsdt_index[1]:dsdt_index[2]],    fmt='k.', capsize=1.5, capthick=1, label='This work (7.8-8.8 GeV)')
plt.errorbar(phi_d_2H_dsdt_minust_center[dsdt_index[2]:dsdt_index[3]],  phi_d_2H_dsdt_results[0, dsdt_index[2]:dsdt_index[3]],  xerr=phi_d_2H_dsdt_minust_width[dsdt_index[2]:dsdt_index[3]],   yerr=phi_d_2H_dsdt_results_statserr[0, dsdt_index[2]:dsdt_index[3]],    fmt='r.', capsize=1.5, capthick=1, label='This work (8.8-10.8 GeV)')
plt.errorbar(phi_d_2H_dsdt_minust_center[dsdt_index[1]:dsdt_index[2]],  phi_d_2H_dsdt_results[0, dsdt_index[1]:dsdt_index[2]],  xerr=phi_d_2H_dsdt_minust_width[dsdt_index[1]:dsdt_index[2]],   yerr=phi_d_2H_dsdt_results_totalerr[dsdt_index[1]:dsdt_index[2]],       fmt='k.')
plt.errorbar(phi_d_2H_dsdt_minust_center[dsdt_index[2]:dsdt_index[3]],  phi_d_2H_dsdt_results[0, dsdt_index[2]:dsdt_index[3]],  xerr=phi_d_2H_dsdt_minust_width[dsdt_index[2]:dsdt_index[3]],   yerr=phi_d_2H_dsdt_results_totalerr[dsdt_index[2]:dsdt_index[3]],       fmt='r.')

plt.errorbar(phi_d_2H_dsdt_theory_minust,                               phi_d_2H_dsdt_theory_10mb,                                                                                              fmt='-', color = 'g', label='E=8.3 GeV $\sigma_{\phi N}=$10 mb (VMD)')
plt.fill_between(phi_d_2H_dsdt_theory_minust,                           phi_d_2H_dsdt_theory_10mb-0.2*(phi_d_2H_dsdt_theory_10mb), phi_d_2H_dsdt_theory_10mb+0.2*(phi_d_2H_dsdt_theory_10mb),   color='g', alpha=0.2)
plt.errorbar(phi_d_2H_dsdt_theory_minust,                               phi_d_2H_dsdt_theory_30mb,                                                                                              fmt='-', color = 'y', label=r'E=8.3 GeV $\sigma_{\phi N}=$30 mb ($b_{\phi N}=8 \rm \ GeV^{-2}$)')
plt.fill_between(phi_d_2H_dsdt_theory_minust,                           phi_d_2H_dsdt_theory_30mb-0.2*(phi_d_2H_dsdt_theory_30mb), phi_d_2H_dsdt_theory_30mb+0.2*(phi_d_2H_dsdt_theory_30mb),   color='y', alpha=0.2)

plt.xlabel(r'$-t\ [GeV^2/c]$', fontsize=16)
plt.ylabel(r'$d\sigma/dt\ [nb/(GeV^2/c)]$', fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.xlim(0, 2)
plt.ylim(1e-1, 1e3)
plt.yscale('log')
plt.legend()
file_pdf.savefig()
plt.close()

# Constraints on sigma_phiN and b_phiN using chi2 fit to ds/dt
print('Performing chi2 analysis on sigma_phiN...')

sphin_list = np.arange(30,40,0.5)
bphin_list = np.arange(10,15,0.5)
normerr_src     = np.sqrt(0.0656**2+0.05**2+0.002**2 + 0.005**2+0.01**2)
normerr_clas    = 0.2
normerr_theory  = 0.2

chi2_array = np.zeros((len(bphin_list), len(sphin_list)))
for i,sphin in enumerate(sphin_list):
    for j,bphin in enumerate(bphin_list):
        print(f'Calculating chi2 for sigma={sphin}, b={bphin}')
        theory_results = np.loadtxt(f'/work/halld2/home/boyu/src_analysis/plot/vm_d/theory/output/nominal/E_8.3_sigma_%.1f_b_%.1f.txt' % (sphin, bphin))
        nuisance_opt = 0
        for nuisance in np.arange(-10, 10.1, 0.01):
            ndf = 0
            this_chi2 = 0
            for k in range(dsdt_index[1], dsdt_index[2]):
                t_val = phi_d_2H_dsdt_minust_center[k]
                if (t_val < 0.24):
                    continue
                data_val = phi_d_2H_dsdt_results[0,k]
                data_err = phi_d_2H_dsdt_results_p2perr[0,k]
                theory_idx = (np.abs(theory_results[:,0] - t_val)).argmin()  # Find the closest theory point
                theory_val = theory_results[theory_idx, 2]
                this_chi2 += ((data_val - (1 + nuisance*np.sqrt(normerr_src**2 + normerr_theory**2))*theory_val)**2)/(data_err**2)
                ndf += 1
            this_chi2 += nuisance**2  # pull term
            ndf -= 2  # 2 parameters: sigma_phiN and b_phiN
            this_chi2 /= ndf
            if this_chi2 < chi2_array[j, i] or chi2_array[j, i] == 0:
                chi2_array[j, i] = this_chi2
                nuisance_opt = nuisance
        if np.abs(nuisance_opt) >= 9.9:
            print(f'Warning: Nuisance parameter at limit for sigma={sphin}, b={bphin}, nuisance={nuisance_opt}')

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

################################################################# OBSERVABLE RESULTS: DECAY ANGULAR DISTRIBUTIONS #############################################################################################################

print('Plotting decay angular distribution results...')

sdme_array = np.zeros((4, len(Wcostheta_index)-1))
error_array = np.zeros((4, len(Wcostheta_index)-1))

# Wcostheta
fig = plt.figure(figsize=(24, 16), dpi=300)
gs = fig.add_gridspec(2, 3, wspace=0, hspace=0)
axs = gs.subplots(sharey=True).flatten()

for i in range(len(Wcostheta_index)-1):
    curve_fit_params, curve_fit_cov = curve_fit(Wcostheta_func, phi_d_2H_Wcostheta_costheta_center[Wcostheta_index[i]:Wcostheta_index[i+1]], phi_d_2H_Wcostheta_results[0, Wcostheta_index[i]:Wcostheta_index[i+1]], sigma=phi_d_2H_Wcostheta_results_totalerr[Wcostheta_index[i]:Wcostheta_index[i+1]], p0=[0.0])
    curve_fit_residuals = phi_d_2H_Wcostheta_results[0, Wcostheta_index[i]:Wcostheta_index[i+1]] - Wcostheta_func(phi_d_2H_Wcostheta_costheta_center[Wcostheta_index[i]:Wcostheta_index[i+1]], curve_fit_params[0])
    reduced_chi2 = np.sum((curve_fit_residuals/phi_d_2H_Wcostheta_results_totalerr[Wcostheta_index[i]:Wcostheta_index[i+1]])**2)/(len(phi_d_2H_Wcostheta_results[0, Wcostheta_index[i]:Wcostheta_index[i+1]])-1)
    sdme_array[0, i] = curve_fit_params[0]
    error_array[0, i] = np.sqrt(curve_fit_cov[0])
    axs[i].errorbar(phi_d_2H_Wcostheta_costheta_center[Wcostheta_index[i]:Wcostheta_index[i+1]], phi_d_2H_Wcostheta_results[0, Wcostheta_index[i]:Wcostheta_index[i+1]], xerr=phi_d_2H_Wcostheta_costheta_width[Wcostheta_index[i]:Wcostheta_index[i+1]], yerr=phi_d_2H_Wcostheta_results_statserr[0, Wcostheta_index[i]:Wcostheta_index[i+1]], fmt='ko', capsize=3, capthick=1, label='This work', markersize=7)
    axs[i].errorbar(phi_d_2H_Wcostheta_costheta_center[Wcostheta_index[i]:Wcostheta_index[i+1]], phi_d_2H_Wcostheta_results[0, Wcostheta_index[i]:Wcostheta_index[i+1]], xerr=phi_d_2H_Wcostheta_costheta_width[Wcostheta_index[i]:Wcostheta_index[i+1]], yerr=phi_d_2H_Wcostheta_results_totalerr[Wcostheta_index[i]:Wcostheta_index[i+1]], fmt='ko')
    axs[i].plot(np.linspace(-1, 1, 100), Wcostheta_func(np.linspace(-1, 1, 100), curve_fit_params[0]), 'b--', label='Fit')
    axs[i].plot(np.linspace(-1, 1, 100), 0.75*(1-np.linspace(-1, 1, 100)**2), 'r--', label='SCHC+NPE')
    axs[i].text(0.05, 0.95, W_label_list[i], transform=axs[i].transAxes, fontsize=16, verticalalignment='top')
    axs[i].text(0.05, 0.85, r'$\rho^0_{00}=%.2f\pm%.2f$' % (curve_fit_params[0], np.sqrt(curve_fit_cov[0])), transform=axs[i].transAxes, fontsize=16, color='b', ha='left', va='top')
    axs[i].text(0.05, 0.75, r'$\chi^2/\mathrm{dof}=%.2f$' % reduced_chi2, transform=axs[i].transAxes, fontsize=16, color='b', ha='left', va='top')
    axs[i].set_xlim(-1, 1)
    axs[i].tick_params(labelsize=16)

axs[0].set_ylabel(r'$W(\cos \vartheta)$', fontsize=16)
axs[3].set_ylabel(r'$W(\cos \vartheta)$', fontsize=16)
axs[0].set_yticks(np.arange(0, 1, 0.2))
axs[0].set_ylim(0, 1)
axs[2].legend()
axs[4].set_xlabel(xlabel_list[1], fontsize=24)
for i in range(3):
    axs[i].set_xticks([])
    axs[i+3].set_xticks(np.arange(-1, 1, 0.25))

file_pdf.savefig()
plt.close()

# Wdecayphi
fig = plt.figure(figsize=(24, 16), dpi=300)
gs = fig.add_gridspec(2, 3, wspace=0, hspace=0)
axs = gs.subplots(sharey=True).flatten()

for i in range(len(Wdecayphi_index)-1):
    curve_fit_params, curve_fit_cov = curve_fit(Wphi_func, phi_d_2H_Wdecayphi_decayphi_center[Wdecayphi_index[i]:Wdecayphi_index[i+1]], phi_d_2H_Wdecayphi_results[0, Wdecayphi_index[i]:Wdecayphi_index[i+1]], sigma=phi_d_2H_Wdecayphi_results_totalerr[Wdecayphi_index[i]:Wdecayphi_index[i+1]], p0=[0.0])
    curve_fit_residuals = phi_d_2H_Wdecayphi_results[0, Wdecayphi_index[i]:Wdecayphi_index[i+1]] - Wphi_func(phi_d_2H_Wdecayphi_decayphi_center[Wdecayphi_index[i]:Wdecayphi_index[i+1]], curve_fit_params[0])
    reduced_chi2 = np.sum((curve_fit_residuals/phi_d_2H_Wdecayphi_results_totalerr[Wdecayphi_index[i]:Wdecayphi_index[i+1]])**2)/(len(phi_d_2H_Wdecayphi_results[0, Wdecayphi_index[i]:Wdecayphi_index[i+1]])-1)
    sdme_array[1, i] = curve_fit_params[0]
    error_array[1, i] = np.sqrt(curve_fit_cov[0])
    axs[i].errorbar(phi_d_2H_Wdecayphi_decayphi_center[Wdecayphi_index[i]:Wdecayphi_index[i+1]], phi_d_2H_Wdecayphi_results[0, Wdecayphi_index[i]:Wdecayphi_index[i+1]], xerr=phi_d_2H_Wdecayphi_decayphi_width[Wdecayphi_index[i]:Wdecayphi_index[i+1]], yerr=phi_d_2H_Wdecayphi_results_totalerr[Wdecayphi_index[i]:Wdecayphi_index[i+1]], fmt='ko')
    axs[i].errorbar(phi_d_2H_Wdecayphi_decayphi_center[Wdecayphi_index[i]:Wdecayphi_index[i+1]], phi_d_2H_Wdecayphi_results[0, Wdecayphi_index[i]:Wdecayphi_index[i+1]], xerr=phi_d_2H_Wdecayphi_decayphi_width[Wdecayphi_index[i]:Wdecayphi_index[i+1]], yerr=phi_d_2H_Wdecayphi_results_statserr[0, Wdecayphi_index[i]:Wdecayphi_index[i+1]], fmt='ko', capsize=3, capthick=1, label='This work', markersize=7)
    axs[i].plot(np.linspace(-180, 180, 360), Wphi_func(np.linspace(-180, 180, 360), curve_fit_params[0]), 'b--', label='Fit')
    axs[i].plot(np.linspace(-180, 180, 360), np.ones(360), 'r--', label='SCHC+NPE')
    axs[i].text(0.05, 0.95, W_label_list[i], transform=axs[i].transAxes, fontsize=16, verticalalignment='top')
    axs[i].text(0.05, 0.85, r'$\rho^0_{1-1}=%.2f\pm%.2f$' % (curve_fit_params[0], np.sqrt(curve_fit_cov[0])), transform=axs[i].transAxes, fontsize=16, color='b', ha='left', va='top')
    axs[i].text(0.05, 0.75, r'$\chi^2/\mathrm{dof}=%.2f$' % reduced_chi2, transform=axs[i].transAxes, fontsize=16, color='b', ha='left', va='top')
    axs[i].set_xlim(-180, 180)
    axs[i].tick_params(labelsize=16)

axs[0].set_ylabel(r'$2\pi W(\varphi)$', fontsize=16)
axs[3].set_ylabel(r'$2\pi W(\varphi)$', fontsize=16)
axs[0].set_yticks(np.arange(0, 2, 0.5))
axs[0].set_ylim(0, 2)
axs[2].legend()
axs[4].set_xlabel(xlabel_list[2], fontsize=24)
for i in range(3):
    axs[i].set_xticks([])
    axs[i+3].set_xticks(np.arange(-180, 180, 45))

file_pdf.savefig()
plt.close()

# Wpolphi
fig = plt.figure(figsize=(24, 16), dpi=300)
gs = fig.add_gridspec(2, 3, wspace=0, hspace=0)
axs = gs.subplots(sharey=True).flatten()

for i in range(len(Wpolphi_index)-1):
    curve_fit_params, curve_fit_cov = curve_fit(Wphi_func, phi_d_2H_Wpolphi_polphi_center[Wpolphi_index[i]:Wpolphi_index[i+1]], phi_d_2H_Wpolphi_results[0, Wpolphi_index[i]:Wpolphi_index[i+1]], sigma=phi_d_2H_Wpolphi_results_totalerr[Wpolphi_index[i]:Wpolphi_index[i+1]], p0=[0.0])
    curve_fit_residuals = phi_d_2H_Wpolphi_results[0, Wpolphi_index[i]:Wpolphi_index[i+1]] - Wphi_func(phi_d_2H_Wpolphi_polphi_center[Wpolphi_index[i]:Wpolphi_index[i+1]], curve_fit_params[0])
    reduced_chi2 = np.sum((curve_fit_residuals/phi_d_2H_Wpolphi_results_totalerr[Wpolphi_index[i]:Wpolphi_index[i+1]])**2)/(len(phi_d_2H_Wpolphi_results[0, Wpolphi_index[i]:Wpolphi_index[i+1]])-1)
    sdme_array[2, i] = curve_fit_params[0]*2/3
    error_array[2, i] = np.sqrt(curve_fit_cov[0])*2/3
    axs[i].errorbar(phi_d_2H_Wpolphi_polphi_center[Wpolphi_index[i]:Wpolphi_index[i+1]], phi_d_2H_Wpolphi_results[0, Wpolphi_index[i]:Wpolphi_index[i+1]], xerr=phi_d_2H_Wpolphi_polphi_width[Wpolphi_index[i]:Wpolphi_index[i+1]], yerr=phi_d_2H_Wpolphi_results_totalerr[Wpolphi_index[i]:Wpolphi_index[i+1]], fmt='ko')
    axs[i].errorbar(phi_d_2H_Wpolphi_polphi_center[Wpolphi_index[i]:Wpolphi_index[i+1]], phi_d_2H_Wpolphi_results[0, Wpolphi_index[i]:Wpolphi_index[i+1]], xerr=phi_d_2H_Wpolphi_polphi_width[Wpolphi_index[i]:Wpolphi_index[i+1]], yerr=phi_d_2H_Wpolphi_results_statserr[0, Wpolphi_index[i]:Wpolphi_index[i+1]], fmt='ko', capsize=3, capthick=1, label='This work', markersize=7)
    axs[i].plot(np.linspace(-180, 180, 360), Wphi_func(np.linspace(-180, 180, 360), curve_fit_params[0]), 'b--', label='Fit')
    axs[i].plot(np.linspace(-180, 180, 360), np.ones(360), 'r--', label='SCHC+NPE')
    axs[i].text(0.05, 0.95, W_label_list[i], transform=axs[i].transAxes, fontsize=16, verticalalignment='top')
    axs[i].text(0.05, 0.85, r'$\rho^0_{1-1}=%.2f\pm%.2f$' % (curve_fit_params[0]*2/3, np.sqrt(curve_fit_cov[0])*2/3), transform=axs[i].transAxes, fontsize=16, color='b', ha='left', va='top')
    axs[i].text(0.05, 0.75, r'$\chi^2/\mathrm{dof}=%.2f$' % reduced_chi2, transform=axs[i].transAxes, fontsize=16, color='b', ha='left', va='top')
    axs[i].set_xlim(-180, 180)
    axs[i].tick_params(labelsize=16)

axs[0].set_ylabel(r'$2\pi W(\Phi)$', fontsize=16)
axs[3].set_ylabel(r'$2\pi W(\Phi)$', fontsize=16)
axs[0].set_yticks(np.arange(0, 2, 0.5))
axs[0].set_ylim(0, 2)
axs[2].legend()
axs[4].set_xlabel(xlabel_list[3], fontsize=24)
for i in range(3):
    axs[i].set_xticks([])
    axs[i+3].set_xticks(np.arange(-180, 180, 45))

file_pdf.savefig()
plt.close()

# Wpsi
fig = plt.figure(figsize=(24, 16), dpi=300)
gs = fig.add_gridspec(2, 3, wspace=0, hspace=0)
axs = gs.subplots(sharey=True).flatten()

polarization_factor = np.array([0.1, 0.1, 0.2, 0.2, 0.036, 0.0036])

for i in range(len(Wpsi_index)-1):
    curve_fit_params, curve_fit_cov = curve_fit(Wphi_func, phi_d_2H_Wpsi_psi_center[Wpsi_index[i]:Wpsi_index[i+1]], phi_d_2H_Wpsi_results[0, Wpsi_index[i]:Wpsi_index[i+1]], sigma=phi_d_2H_Wpsi_results_totalerr[Wpsi_index[i]:Wpsi_index[i+1]], p0=[0.0])
    curve_fit_residuals = phi_d_2H_Wpsi_results[0, Wpsi_index[i]:Wpsi_index[i+1]] - Wphi_func(phi_d_2H_Wpsi_psi_center[Wpsi_index[i]:Wpsi_index[i+1]], curve_fit_params[0])
    reduced_chi2 = np.sum((curve_fit_residuals/phi_d_2H_Wpsi_results_totalerr[Wpsi_index[i]:Wpsi_index[i+1]])**2)/(len(phi_d_2H_Wpsi_results[0, Wpsi_index[i]:Wpsi_index[i+1]])-1)
    sdme_array[3, i] = -curve_fit_params[0]/polarization_factor[i]
    error_array[3, i] = np.sqrt(curve_fit_cov[0])/polarization_factor[i]
    axs[i].errorbar(phi_d_2H_Wpsi_psi_center[Wpsi_index[i]:Wpsi_index[i+1]], phi_d_2H_Wpsi_results[0, Wpsi_index[i]:Wpsi_index[i+1]], xerr=phi_d_2H_Wpsi_psi_width[Wpsi_index[i]:Wpsi_index[i+1]], yerr=phi_d_2H_Wpsi_results_totalerr[Wpsi_index[i]:Wpsi_index[i+1]], fmt='ko')
    axs[i].errorbar(phi_d_2H_Wpsi_psi_center[Wpsi_index[i]:Wpsi_index[i+1]], phi_d_2H_Wpsi_results[0, Wpsi_index[i]:Wpsi_index[i+1]], xerr=phi_d_2H_Wpsi_psi_width[Wpsi_index[i]:Wpsi_index[i+1]], yerr=phi_d_2H_Wpsi_results_statserr[0, Wpsi_index[i]:Wpsi_index[i+1]], fmt='ko', capsize=3, capthick=1, label='This work', markersize=7)
    axs[i].plot(np.linspace(-180, 180, 360), Wphi_func(np.linspace(-180, 180, 360), curve_fit_params[0]), 'b--', label='Fit')
    axs[i].plot(np.linspace(-180, 180, 360), np.ones(360)+polarization_factor[i]*np.cos(2*np.linspace(-180, 180, 360)/rad_to_deg), 'r--', label='SCHC+NPE')
    axs[i].text(0.05, 0.95, W_label_list[i], transform=axs[i].transAxes, fontsize=16, verticalalignment='top')
    axs[i].text(0.05, 0.85, r'$\rho^0_{1-1}=%.2f\pm%.2f$' % (-curve_fit_params[0]/polarization_factor[i], np.sqrt(curve_fit_cov[0])/polarization_factor[i]), transform=axs[i].transAxes, fontsize=16, color='b', ha='left', va='top')
    axs[i].text(0.05, 0.75, r'$\chi^2/\mathrm{dof}=%.2f$' % reduced_chi2, transform=axs[i].transAxes, fontsize=16, color='b', ha='left', va='top')
    axs[i].set_xlim(-180, 180)
    axs[i].tick_params(labelsize=16)

axs[0].set_ylabel(r'$2\pi W(\psi)$', fontsize=16)
axs[3].set_ylabel(r'$2\pi W(\psi)$', fontsize=16)
axs[0].set_yticks(np.arange(0, 2, 0.5))
axs[0].set_ylim(0, 2)
axs[2].legend()
axs[4].set_xlabel(xlabel_list[4], fontsize=24)
for i in range(3):
    axs[i].set_xticks([])
    axs[i+3].set_xticks(np.arange(-180, 180, 45))

file_pdf.savefig()
plt.close()

# SDMEs
fig = plt.figure(figsize=(16, 16), dpi=300)
gs = fig.add_gridspec(2, 2, wspace=0.2, hspace=0.2)
axs = gs.subplots().flatten()

xpoint_list = np.array([6.8739, 6.8739+0.1, 8.3011, 8.3011+0.1, 9.6774, 9.6774+0.1])

for i in range(3):
    axs[i].set_ylim(-0.2, 0.2)
    axs[i].plot(np.linspace(6, 11, 100), np.zeros(100), 'r--', label='SCHC+NPE')
axs[3].set_ylim(0.0, 1.0)
axs[3].plot(np.linspace(6, 11, 100), 0.5*np.ones(100), 'r--', label='SCHC+NPE')
axs[0].set_ylabel(r'$\rho^0_{00}$', fontsize=16)
axs[1].set_ylabel(r'$\rho^0_{1-1}$', fontsize=16)
axs[2].set_ylabel(r'$\rho^1_{00}$', fontsize=16)
axs[3].set_ylabel(r'$\rho^1_{1-1}$', fontsize=16)

for i in range(4):
    for j in range(3):
        axs[i].errorbar(xpoint_list[2*j], sdme_array[i, 2*j], yerr=error_array[i, 2*j], fmt='ko', capsize=3, capthick=1, label=r"$-t < 0.6 \rm{ GeV}^2$", markersize=7)
        axs[i].errorbar(xpoint_list[2*j+1], sdme_array[i, 2*j+1], yerr=error_array[i, 2*j+1], fmt='bo', capsize=3, capthick=1, label=r"$-t > 0.6 \rm{ GeV}^2$", markersize=7)
    axs[i].set_xlim(6, 11)
    axs[i].tick_params(labelsize=16)
    axs[i].set_xlabel("Photon energy (GeV)", fontsize=16)
    axs[i].legend()
    legend_without_duplicate_labels(axs[i])

#     curve_fit_params, curve_fit_cov = curve_fit(Wphi_func, phi_d_2H_Wpsi_psi_center[Wpsi_index[i]:Wpsi_index[i+1]], phi_d_2H_Wpsi_results[0, Wpsi_index[i]:Wpsi_index[i+1]], sigma=phi_d_2H_Wpsi_results_totalerr[Wpsi_index[i]:Wpsi_index[i+1]], p0=[0.0])
#     curve_fit_residuals = phi_d_2H_Wpsi_results[0, Wpsi_index[i]:Wpsi_index[i+1]] - Wphi_func(phi_d_2H_Wpsi_psi_center[Wpsi_index[i]:Wpsi_index[i+1]], curve_fit_params[0])
#     reduced_chi2 = np.sum((curve_fit_residuals/phi_d_2H_Wpsi_results_totalerr[Wpsi_index[i]:Wpsi_index[i+1]])**2)/(len(phi_d_2H_Wpsi_results[0, Wpsi_index[i]:Wpsi_index[i+1]])-1)
#     sdme_array[3, i] = curve_fit_params[0]/polarization_factor[i]
#     error_array[3, i] = np.sqrt(curve_fit_cov[0])/polarization_factor[i]
#     axs[i].errorbar(phi_d_2H_Wpsi_psi_center[Wpsi_index[i]:Wpsi_index[i+1]], phi_d_2H_Wpsi_results[0, Wpsi_index[i]:Wpsi_index[i+1]], xerr=phi_d_2H_Wpsi_psi_width[Wpsi_index[i]:Wpsi_index[i+1]], yerr=phi_d_2H_Wpsi_results_totalerr[Wpsi_index[i]:Wpsi_index[i+1]], fmt='ko')
#     axs[i].errorbar(phi_d_2H_Wpsi_psi_center[Wpsi_index[i]:Wpsi_index[i+1]], phi_d_2H_Wpsi_results[0, Wpsi_index[i]:Wpsi_index[i+1]], xerr=phi_d_2H_Wpsi_psi_width[Wpsi_index[i]:Wpsi_index[i+1]], yerr=phi_d_2H_Wpsi_results_statserr[0, Wpsi_index[i]:Wpsi_index[i+1]], fmt='ko', capsize=3, capthick=1, label='This work', markersize=7)
#     axs[i].plot(np.linspace(-180, 180, 360), Wphi_func(np.linspace(-180, 180, 360), curve_fit_params[0]), 'b--', label='Fit')
#     axs[i].plot(np.linspace(-180, 180, 360), np.ones(360)+polarization_factor[i]*np.cos(2*np.linspace(-180, 180, 360)/rad_to_deg), 'r--', label='SCHC+NPE')
#     axs[i].text(0.05, 0.95, W_label_list[i], transform=axs[i].transAxes, fontsize=16, verticalalignment='top')
#     axs[i].text(0.05, 0.85, r'$\rho^0_{1-1}=%.2f\pm%.2f$' % (-curve_fit_params[0]/polarization_factor[i], np.sqrt(curve_fit_cov[0])/polarization_factor[i]), transform=axs[i].transAxes, fontsize=16, color='b', ha='left', va='top')
#     axs[i].text(0.05, 0.75, r'$\chi^2/\mathrm{dof}=%.2f$' % reduced_chi2, transform=axs[i].transAxes, fontsize=16, color='b', ha='left', va='top')
#     axs[i].set_xlim(-180, 180)
#     axs[i].tick_params(labelsize=16)

# axs[0].set_ylabel(r'$2\pi W(\psi)$', fontsize=16)
# axs[3].set_ylabel(r'$2\pi W(\psi)$', fontsize=16)
# axs[0].set_yticks(np.arange(0, 2, 0.5))
# axs[0].set_ylim(0, 2)
# axs[2].legend()
# axs[4].set_xlabel(xlabel_list[4], fontsize=24)
# for i in range(3):
#     axs[i].set_xticks([])
#     axs[i+3].set_xticks(np.arange(-180, 180, 45))

file_pdf.savefig()
plt.close()

################################################################# OUTPUT RESULT FOR LaTeX #############################################################################################################

print('Writing results to text file...')

file_txt = open("output/table_numerical.txt", "w")

def format_scientific(number, precision=2):
    # Format the number in scientific notation (e.g., 'e' or 'E')
    formatted = f"{number:.{precision}E}"
    # Remove the '+', but keep the '-' for negative exponents
    return formatted.replace('e+0', 'e').replace('E+0', 'E')


file_txt.write("dsdt\n")
for i in range(phi_d_2H_dsdt_results[0].shape[0]):
    file_txt.write("        ")
    file_txt.write(f"{phi_d_2H_dsdt_minust_low[i]:5.3f} & {phi_d_2H_dsdt_minust_high[i]:5.3f} & ")
    file_txt.write(f"{phi_d_2H_dsdt_minust_center[i]:5.4f} & {phi_d_2H_dsdt_minust_width[i]:5.4f} &    ")
    file_txt.write(f"{phi_d_2H_dsdt_yield_data[0, i]:6.2f} $\pm$ {phi_d_2H_dsdt_yield_data_statserr[0, i]:5.2f} &    ")
    file_txt.write(f"{format_scientific(phi_d_2H_dsdt_yield_sim[0, i])} $\pm$ {format_scientific(phi_d_2H_dsdt_yield_sim_statserr[0, i])} &    ")
    file_txt.write(f"{format_scientific(phi_d_2H_dsdt_yield_tagged[0, i])} $\pm$ {format_scientific(phi_d_2H_dsdt_yield_tagged_statserr[0, i])} &    ")
    file_txt.write(f"{phi_d_2H_dsdt_results[0, i]:6.2f} $\pm$ {phi_d_2H_dsdt_results_statserr[0, i]:5.2f} $\pm$ {phi_d_2H_dsdt_results_systerr[i]:5.2f}")
    file_txt.write(" \\\\\n")
    if (i+1) in dsdt_index:
        file_txt.write("\n")

file_txt.write("Wcostheta\n")
for i in range(phi_d_2H_Wcostheta_results[0].shape[0]):
    file_txt.write("        & ")
    file_txt.write(f"{phi_d_2H_Wcostheta_costheta_low[i]:5.1f} & {phi_d_2H_Wcostheta_costheta_high[i]:5.1f} & ")
    file_txt.write(f"{phi_d_2H_Wcostheta_costheta_center[i]:6.3f} & {phi_d_2H_Wcostheta_costheta_width[i]:5.3f} &    ")
    file_txt.write(f"{phi_d_2H_Wcostheta_yield_data[0, i]:6.2f} $\pm$ {phi_d_2H_Wcostheta_yield_data_statserr[0, i]:5.2f} &    ")
    file_txt.write(f"{format_scientific(phi_d_2H_Wcostheta_yield_sim[0, i])} $\pm$ {format_scientific(phi_d_2H_Wcostheta_yield_sim_statserr[0, i])} &    ")
    file_txt.write(f"{format_scientific(phi_d_2H_Wcostheta_yield_tagged[0, i])} $\pm$ {format_scientific(phi_d_2H_Wcostheta_yield_tagged_statserr[0, i])} &    ")
    file_txt.write(f"{phi_d_2H_Wcostheta_results[0, i]:6.2f} $\pm$ {phi_d_2H_Wcostheta_results_statserr[0, i]:4.2f} $\pm$ {phi_d_2H_Wcostheta_results_systerr[i]:4.2f}")
    file_txt.write(" \\\\\n")
    if (i+1) in Wcostheta_index:
        file_txt.write("\n")

file_txt.write("Wdecayphi\n")
for i in range(phi_d_2H_Wdecayphi_results[0].shape[0]):
    file_txt.write("        & ")
    file_txt.write(f"{phi_d_2H_Wdecayphi_decayphi_low[i]:4.0f} & {phi_d_2H_Wdecayphi_decayphi_high[i]:4.0f} & ")
    file_txt.write(f"{phi_d_2H_Wdecayphi_decayphi_center[i]:7.2f} & {phi_d_2H_Wdecayphi_decayphi_width[i]:5.2f} &    ")
    file_txt.write(f"{phi_d_2H_Wdecayphi_yield_data[0, i]:6.2f} $\pm$ {phi_d_2H_Wdecayphi_yield_data_statserr[0, i]:5.2f} &    ")
    file_txt.write(f"{format_scientific(phi_d_2H_Wdecayphi_yield_sim[0, i])} $\pm$ {format_scientific(phi_d_2H_Wdecayphi_yield_sim_statserr[0, i])} &    ")
    file_txt.write(f"{format_scientific(phi_d_2H_Wdecayphi_yield_tagged[0, i])} $\pm$ {format_scientific(phi_d_2H_Wdecayphi_yield_tagged_statserr[0, i])} &    ")
    file_txt.write(f"{phi_d_2H_Wdecayphi_results[0, i]:6.2f} $\pm$ {phi_d_2H_Wdecayphi_results_statserr[0, i]:4.2f} $\pm$ {phi_d_2H_Wdecayphi_results_systerr[i]:4.2f}")
    file_txt.write(" \\\\\n")
    if (i+1) in Wdecayphi_index:
        file_txt.write("\n")

file_txt.write("Wpolphi\n")
for i in range(phi_d_2H_Wpolphi_results[0].shape[0]):
    file_txt.write("        & ")
    file_txt.write(f"{phi_d_2H_Wpolphi_polphi_low[i]:4.0f} & {phi_d_2H_Wpolphi_polphi_high[i]:4.0f} & ")
    file_txt.write(f"{phi_d_2H_Wpolphi_polphi_center[i]:7.2f} & {phi_d_2H_Wpolphi_polphi_width[i]:5.2f} &    ")
    file_txt.write(f"{phi_d_2H_Wpolphi_yield_data[0, i]:6.2f} $\pm$ {phi_d_2H_Wpolphi_yield_data_statserr[0, i]:5.2f} &    ")
    file_txt.write(f"{format_scientific(phi_d_2H_Wpolphi_yield_sim[0, i])} $\pm$ {format_scientific(phi_d_2H_Wpolphi_yield_sim_statserr[0, i])} &    ")
    file_txt.write(f"{format_scientific(phi_d_2H_Wpolphi_yield_tagged[0, i])} $\pm$ {format_scientific(phi_d_2H_Wpolphi_yield_tagged_statserr[0, i])} &    ")
    file_txt.write(f"{phi_d_2H_Wpolphi_results[0, i]:6.2f} $\pm$ {phi_d_2H_Wpolphi_results_statserr[0, i]:4.2f} $\pm$ {phi_d_2H_Wpolphi_results_systerr[i]:4.2f}")
    file_txt.write(" \\\\\n")
    if (i+1) in Wpolphi_index:
        file_txt.write("\n")

file_txt.write("Wpsi\n")
for i in range(phi_d_2H_Wpsi_results[0].shape[0]):
    file_txt.write("        & ")
    file_txt.write(f"{phi_d_2H_Wpsi_psi_low[i]:4.0f} & {phi_d_2H_Wpsi_psi_high[i]:4.0f} & ")
    file_txt.write(f"{phi_d_2H_Wpsi_psi_center[i]:7.2f} & {phi_d_2H_Wpsi_psi_width[i]:5.2f} &    ")
    file_txt.write(f"{phi_d_2H_Wpsi_yield_data[0, i]:6.2f} $\pm$ {phi_d_2H_Wpsi_yield_data_statserr[0, i]:5.2f} &    ")
    file_txt.write(f"{format_scientific(phi_d_2H_Wpsi_yield_sim[0, i])} $\pm$ {format_scientific(phi_d_2H_Wpsi_yield_sim_statserr[0, i])} &    ")
    file_txt.write(f"{format_scientific(phi_d_2H_Wpsi_yield_tagged[0, i])} $\pm$ {format_scientific(phi_d_2H_Wpsi_yield_tagged_statserr[0, i])} &    ")
    file_txt.write(f"{phi_d_2H_Wpsi_results[0, i]:6.2f} $\pm$ {phi_d_2H_Wpsi_results_statserr[0, i]:4.2f} $\pm$ {phi_d_2H_Wpsi_results_systerr[i]:4.2f}")
    file_txt.write(" \\\\\n")
    if (i+1) in Wpsi_index:
        file_txt.write("\n")

file_txt.close()





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