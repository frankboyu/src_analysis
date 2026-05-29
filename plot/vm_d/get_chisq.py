import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages

rad_to_deg = 180/np.pi
# file_pdf = PdfPages("/work/halld2/home/boyu/src_analysis/plot/vm_d/output/plots_vm_d.pdf")

def legend_without_duplicate_labels(ax):
    handles, labels = ax.get_legend_handles_labels()
    unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
    ax.legend(*zip(*unique))

################################################################# DATA POINTS ###################################################################################

# Normalization uncertainties
normerr_theory      = 0.2
normerr_src         = np.sqrt(0.0656**2+0.05**2+0.002**2 + 0.005**2+0.01**2)
normerr_src_total   = np.sqrt(normerr_src**2 + normerr_theory**2)
normerr_clas        = np.sqrt(0.005**2 + 0.10**2 + 0.005**2 + 0.003**2 + 0.014**2)
normerr_clas_total  = np.sqrt(normerr_clas**2 + normerr_theory**2)

# Load ds/dt results
phi_d_2H_dsdt_energy_center     = np.loadtxt('output/table_vm_d_dsdt.txt')[:,0]
phi_d_2H_dsdt_minust_center     = np.loadtxt('output/table_vm_d_dsdt.txt')[:,1]
phi_d_2H_dsdt_results           = np.loadtxt('output/table_vm_d_dsdt.txt')[:,2]
phi_d_2H_dsdt_results_statserr  = np.loadtxt('output/table_vm_d_dsdt.txt')[:,3]
phi_d_2H_dsdt_results_p2perr    = np.loadtxt('output/table_vm_d_dsdt.txt')[:,4]
phi_d_2H_dsdt_results_normerr   = np.loadtxt('output/table_vm_d_dsdt.txt')[:,5]
phi_d_2H_dsdt_results_systerr   = np.loadtxt('output/table_vm_d_dsdt.txt')[:,6]

# Identify indices for different energy bins based on gaps in energy_center
dsdt_index = []
for i in range(len(phi_d_2H_dsdt_results)):
    if (i == 0):
        dsdt_index.append(i)
    elif (i == len(phi_d_2H_dsdt_results) - 1):
        dsdt_index.append(i+1)
    else:
        if (phi_d_2H_dsdt_energy_center[i]-phi_d_2H_dsdt_energy_center[i-1] > 1.0):
            dsdt_index.append(i)

# CLAS results
phi_d_2H_dsdt_clas_minust_center        = np.array([0.360,  0.385,  0.410,  0.435,  0.474,  0.524,  0.574,  0.646,  0.746,  0.888,  1.091,  1.292,  1.637])
phi_d_2H_dsdt_clas_results_16           = np.array([10.21,  8.85,   7.32,   6.16,   4.73,   3.52,   2.66,   2.17,   1.40,   0.94,   0.57,   0.28,   0.19])
phi_d_2H_dsdt_clas_results_16_statserr  = np.array([0.82,   0.75,   0.59,   0.55,   0.34,   0.28,   0.24,   0.15,   0.12,   0.07,   0.06,   0.05,   0.02])
phi_d_2H_dsdt_clas_results_16_p2perr    = np.array([13.0,   7.2,    7.9,    8.2,    7.3,    10.1,   9.6,    6.2,    5.1,    6.4,    5.9,    9.3,    9.9])/100.0 * phi_d_2H_dsdt_clas_results_16
phi_d_2H_dsdt_clas_results_26           = np.array([8.63,   6.80,   4.57,   5.76,   3.99,   3.59,   2.11,   1.83,   1.32,   0.96,   0.57,   0.36,   0.15])
phi_d_2H_dsdt_clas_results_26_statserr  = np.array([0.80,   0.69,   0.53,   0.56,   0.33,   0.29,   0.22,   0.14,   0.12,   0.07,   0.05,   0.04,   0.02])
phi_d_2H_dsdt_clas_results_26_p2perr    = np.array([6.4,    11.9,   12.3,   4.8,    9.0,    11.1,   8.1,    8.5,    10.7,   5.4,    4.2,    6.9,    4.3])/100.0 * phi_d_2H_dsdt_clas_results_26

# Array to store the results of chi2 analysis
sphin_list      = np.arange(30, 40, 1)
bphin_list      = np.arange(10, 15, 0.25)
chi2_array      = np.zeros((8, len(sphin_list), len(bphin_list)))
ndf_array       = np.zeros((8, len(sphin_list), len(bphin_list)))
nuisance_array  = np.zeros((8, len(sphin_list), len(bphin_list)))

# ################################################################# SRC-CT 6.9 GEV ###################################################################################

# print('Calculating chi2 for SRC-CT 6.9 GeV data...')
# index = 0

# for i,sphin in enumerate(sphin_list):
#     print(f'sigma={sphin}')
#     for j,bphin in enumerate(bphin_list):
#         theory_results = np.loadtxt(f'/work/halld2/home/boyu/src_analysis/plot/vm_d/theory/output/nominal/E_6.9_sigma_%.0f_b_%.2f.txt' % (sphin, bphin))
#         for this_nuisance in np.arange(-10, 10.1, 0.01):
#             this_ndf = 0
#             this_chi2 = 0
#             for k in range(dsdt_index[index], dsdt_index[index+1]):
#                 t_val = phi_d_2H_dsdt_minust_center[k]
#                 data_val = phi_d_2H_dsdt_results[k]
#                 data_err = np.sqrt(phi_d_2H_dsdt_results_p2perr[k]**2 + phi_d_2H_dsdt_results_statserr[k]**2)
#                 theory_idx = (np.abs(theory_results[:,0] - t_val)).argmin()  # Find the closest theory point
#                 theory_val = theory_results[theory_idx, 2]
#                 this_chi2 += ((data_val - (1 + this_nuisance*normerr_src_total)*theory_val)**2)/(data_err**2)
#                 this_ndf += 1
#             this_chi2 += this_nuisance**2  # pull term
#             this_ndf -= 2  # 2 parameters: sigma_phiN and b_phiN
#             if this_chi2 < chi2_array[index, i, j] or chi2_array[index, i, j] == 0:
#                 chi2_array[index, i, j]     = this_chi2
#                 ndf_array[index, i, j]      = this_ndf
#                 nuisance_array[index, i, j] = this_nuisance

# ################################################################# SRC-CT 8.3 GEV ###################################################################################

# print('Calculating chi2 for SRC-CT 8.3 GeV data...')
# index = 1

# for i,sphin in enumerate(sphin_list):
#     print(f'sigma={sphin}')
#     for j,bphin in enumerate(bphin_list):
#         theory_results = np.loadtxt(f'/work/halld2/home/boyu/src_analysis/plot/vm_d/theory/output/nominal/E_8.3_sigma_%.0f_b_%.2f.txt' % (sphin, bphin))
#         for this_nuisance in np.arange(-10, 10.1, 0.01):
#             this_ndf = 0
#             this_chi2 = 0
#             for k in range(dsdt_index[index], dsdt_index[index+1]):
#                 t_val = phi_d_2H_dsdt_minust_center[k]
#                 data_val = phi_d_2H_dsdt_results[k]
#                 data_err = np.sqrt(phi_d_2H_dsdt_results_p2perr[k]**2 + phi_d_2H_dsdt_results_statserr[k]**2)
#                 theory_idx = (np.abs(theory_results[:,0] - t_val)).argmin()  # Find the closest theory point
#                 theory_val = theory_results[theory_idx, 2]
#                 this_chi2 += ((data_val - (1 + this_nuisance*normerr_src_total)*theory_val)**2)/(data_err**2)
#                 this_ndf += 1
#             this_chi2 += this_nuisance**2  # pull term
#             this_ndf -= 2  # 2 parameters: sigma_phiN and b_phiN
#             if this_chi2 < chi2_array[index, i, j] or chi2_array[index, i, j] == 0:
#                 chi2_array[index, i, j]     = this_chi2
#                 ndf_array[index, i, j]      = this_ndf
#                 nuisance_array[index, i, j] = this_nuisance

# ################################################################# SRC-CT 9.7 GEV ###################################################################################

# print('Calculating chi2 for SRC-CT 9.7 GeV data...')
# index = 2

# for i,sphin in enumerate(sphin_list):
#     print(f'sigma={sphin}')
#     for j,bphin in enumerate(bphin_list):
#         theory_results = np.loadtxt(f'/work/halld2/home/boyu/src_analysis/plot/vm_d/theory/output/nominal/E_9.7_sigma_%.0f_b_%.2f.txt' % (sphin, bphin))
#         for this_nuisance in np.arange(-10, 10.1, 0.01):
#             this_ndf = 0
#             this_chi2 = 0
#             for k in range(dsdt_index[index], dsdt_index[index+1]):
#                 t_val = phi_d_2H_dsdt_minust_center[k]
#                 data_val = phi_d_2H_dsdt_results[k]
#                 data_err = np.sqrt(phi_d_2H_dsdt_results_p2perr[k]**2 + phi_d_2H_dsdt_results_statserr[k]**2)
#                 theory_idx = (np.abs(theory_results[:,0] - t_val)).argmin()  # Find the closest theory point
#                 theory_val = theory_results[theory_idx, 2]
#                 this_chi2 += ((data_val - (1 + this_nuisance*normerr_src_total)*theory_val)**2)/(data_err**2)
#                 this_ndf += 1
#             this_chi2 += this_nuisance**2  # pull term
#             this_ndf -= 2  # 2 parameters: sigma_phiN and b_phiN
#             if this_chi2 < chi2_array[index, i, j] or chi2_array[index, i, j] == 0:
#                 chi2_array[index, i, j]     = this_chi2
#                 ndf_array[index, i, j]      = this_ndf
#                 nuisance_array[index, i, j] = this_nuisance

# ################################################################# SRC-CT COMBINED ###################################################################################

# print('Calculating chi2 for SRC-CT data combined...')
# index = 3

# for i,sphin in enumerate(sphin_list):
#     print(f'sigma={sphin}')
#     for j,bphin in enumerate(bphin_list):
#         theory_results = []
#         theory_results.append(np.loadtxt(f'/work/halld2/home/boyu/src_analysis/plot/vm_d/theory/output/nominal/E_6.9_sigma_%.0f_b_%.2f.txt' % (sphin, bphin)))
#         theory_results.append(np.loadtxt(f'/work/halld2/home/boyu/src_analysis/plot/vm_d/theory/output/nominal/E_8.3_sigma_%.0f_b_%.2f.txt' % (sphin, bphin)))
#         theory_results.append(np.loadtxt(f'/work/halld2/home/boyu/src_analysis/plot/vm_d/theory/output/nominal/E_9.7_sigma_%.0f_b_%.2f.txt' % (sphin, bphin)))
#         for this_nuisance in np.arange(-10, 10.1, 0.01):
#             this_ndf = 0
#             this_chi2 = 0
#             for l in range(3):
#                 for k in range(dsdt_index[l], dsdt_index[l+1]):
#                     t_val = phi_d_2H_dsdt_minust_center[k]
#                     data_val = phi_d_2H_dsdt_results[k]
#                     data_err = np.sqrt(phi_d_2H_dsdt_results_p2perr[k]**2 + phi_d_2H_dsdt_results_statserr[k]**2)
#                     theory_idx = (np.abs(theory_results[l][:,0] - t_val)).argmin()  # Find the closest theory point
#                     theory_val = theory_results[l][theory_idx, 2]
#                     this_chi2 += ((data_val - (1 + this_nuisance*normerr_src_total)*theory_val)**2)/(data_err**2)
#                     this_ndf += 1
#             this_chi2 += this_nuisance**2  # pull term
#             this_ndf -= 2  # 2 parameters: sigma_phiN and b_phiN
#             if this_chi2 < chi2_array[index, i, j] or chi2_array[index, i, j] == 0:
#                 chi2_array[index, i, j]     = this_chi2
#                 ndf_array[index, i, j]      = this_ndf
#                 nuisance_array[index, i, j] = this_nuisance

################################################################# CLAS 2.1 GEV ###################################################################################

print('Calculating chi2 for CLAS 2.1 GeV data...')
index = 4

for i,sphin in enumerate(sphin_list):
    print(f'sigma={sphin}')
    for j,bphin in enumerate(bphin_list):
        theory_results = np.loadtxt(f'/work/halld2/home/boyu/src_analysis/plot/vm_d/theory/output/nominal/E_2.1_sigma_%.0f_b_%.2f.txt' % (sphin, bphin))
        for this_nuisance in np.arange(-10, 10.1, 0.01):
            this_ndf    = 0
            this_chi2   = 0
            for k in range(len(phi_d_2H_dsdt_clas_minust_center)):
                t_val       = phi_d_2H_dsdt_clas_minust_center[k]
                data_val    = phi_d_2H_dsdt_clas_results_16[k]
                data_err    = np.sqrt(phi_d_2H_dsdt_clas_results_16_p2perr[k]**2 + phi_d_2H_dsdt_clas_results_16_statserr[k]**2)
                theory_idx  = (np.abs(theory_results[:,0] - t_val)).argmin()  # Find the closest theory point
                theory_val  = theory_results[theory_idx, 2]
                this_chi2   += ((data_val - (1 + this_nuisance*normerr_clas_total)*theory_val)**2)/(data_err**2)
                this_ndf    += 1
            this_chi2   += this_nuisance**2     # pull term
            this_ndf    -= 2                    # 2 parameters: sigma_phiN and b_phiN
            if this_chi2 < chi2_array[index, i, j] or chi2_array[index, i, j] == 0:
                chi2_array[index, i, j]     = this_chi2
                ndf_array[index, i, j]      = this_ndf
                nuisance_array[index, i, j] = this_nuisance

################################################################# CLAS 3.1 GEV ###################################################################################

print('Calculating chi2 for CLAS 3.1 GeV data...')
index = 5

for i,sphin in enumerate(sphin_list):
    print(f'sigma={sphin}')
    for j,bphin in enumerate(bphin_list):
        theory_results = np.loadtxt(f'/work/halld2/home/boyu/src_analysis/plot/vm_d/theory/output/nominal/E_3.1_sigma_%.0f_b_%.2f.txt' % (sphin, bphin))
        for this_nuisance in np.arange(-10, 10.1, 0.01):
            this_ndf = 0
            this_chi2 = 0
            for k in range(len(phi_d_2H_dsdt_clas_minust_center)):
                t_val       = phi_d_2H_dsdt_clas_minust_center[k]
                data_val    = phi_d_2H_dsdt_clas_results_16[k]
                data_err    = np.sqrt(phi_d_2H_dsdt_clas_results_16_p2perr[k]**2 + phi_d_2H_dsdt_clas_results_16_statserr[k]**2)
                theory_idx  = (np.abs(theory_results[:,0] - t_val)).argmin()  # Find the closest theory point
                theory_val  = theory_results[theory_idx, 2]
                this_chi2   += ((data_val - (1 + this_nuisance*normerr_clas_total)*theory_val)**2)/(data_err**2)
                this_ndf    += 1
            this_chi2   += this_nuisance**2     # pull term
            this_ndf    -= 2                    # 2 parameters: sigma_phiN and b_phiN
            if this_chi2 < chi2_array[index, i, j] or chi2_array[index, i, j] == 0:
                chi2_array[index, i, j]     = this_chi2
                ndf_array[index, i, j]      = this_ndf
                nuisance_array[index, i, j] = this_nuisance

################################################################# CLAS COMBINED ###################################################################################

print('Calculating chi2 for CLAS data combined...')
index = 6

for i,sphin in enumerate(sphin_list):
    print(f'sigma={sphin}')
    for j,bphin in enumerate(bphin_list):
        theory_results = []
        theory_results.append(np.loadtxt(f'/work/halld2/home/boyu/src_analysis/plot/vm_d/theory/output/nominal/E_2.1_sigma_%.0f_b_%.2f.txt' % (sphin, bphin)))
        theory_results.append(np.loadtxt(f'/work/halld2/home/boyu/src_analysis/plot/vm_d/theory/output/nominal/E_3.1_sigma_%.0f_b_%.2f.txt' % (sphin, bphin)))
        for this_nuisance in np.arange(-10, 10.1, 0.01):
            this_ndf = 0
            this_chi2 = 0
            for k in range(len(phi_d_2H_dsdt_clas_minust_center)):
                t_val = phi_d_2H_dsdt_clas_minust_center[k]
                data_val = phi_d_2H_dsdt_clas_results_16[k]
                data_err = np.sqrt(phi_d_2H_dsdt_clas_results_16_p2perr[k]**2 + phi_d_2H_dsdt_clas_results_16_statserr[k]**2)
                theory_idx = (np.abs(theory_results[0][:,0] - t_val)).argmin()  # Find the closest theory point
                theory_val = theory_results[0][theory_idx, 2]
                this_chi2 += ((data_val - (1 + this_nuisance*normerr_clas_total)*theory_val)**2)/(data_err**2)
                this_ndf += 1
            for k in range(len(phi_d_2H_dsdt_clas_minust_center)):
                t_val       = phi_d_2H_dsdt_clas_minust_center[k]
                data_val    = phi_d_2H_dsdt_clas_results_26[k]
                data_err    = np.sqrt(phi_d_2H_dsdt_clas_results_26_p2perr[k]**2 + phi_d_2H_dsdt_clas_results_26_statserr[k]**2)
                theory_idx  = (np.abs(theory_results[1][:,0] - t_val)).argmin()  # Find the closest theory point
                theory_val  = theory_results[1][theory_idx, 2]
                this_chi2   += ((data_val - (1 + this_nuisance*normerr_clas_total)*theory_val)**2)/(data_err**2)
                this_ndf    += 1
            this_chi2   += this_nuisance**2     # pull term
            this_ndf    -= 2                    # 2 parameters: sigma_phiN and b_phiN
            if this_chi2 < chi2_array[index, i, j] or chi2_array[index, i, j] == 0:
                chi2_array[index, i, j]     = this_chi2
                ndf_array[index, i, j]      = this_ndf
                nuisance_array[index, i, j] = this_nuisance

################################################################# SRC-CT AND CLAS COMBINED ###################################################################################

################################################################# SAVE RESULTS ###################################################################################

file_txt = open('output/table_vm_d_chisq.txt', 'w')

for i,sphin in enumerate(sphin_list):
    for j,bphin in enumerate(bphin_list):
        file_txt.write(f"{sphin:5.1f} \t {bphin:5.2f} \t \t")
        for k in range(8):
            file_txt.write(f"\t {chi2_array[k, i, j]:8.4f} \t {ndf_array[k, i, j]:4.0f} \t {nuisance_array[k, i, j]:7.2f}")
            if (k != 7):
                file_txt.write("\t \t")
        file_txt.write("\n")
