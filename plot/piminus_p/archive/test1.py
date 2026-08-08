# Import packages
import numpy as np
import ROOT as root
import matplotlib
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages

def create_yield_hist():
    yield_hist = []
    for i in range(len(edges_energy)-1):
        yield_hist_energy = []
        for j in range(len(edges_theta)-1):
            yield_hist_theta = root.TH1F("", "", 50, 0.5, 1.3)
            yield_hist_energy.append(yield_hist_theta)
        yield_hist.append(yield_hist_energy)
    return yield_hist

def hist_points_1D(hist):
    Nbins = hist.GetNbinsX()
    x, y, x_err, y_err = np.empty(Nbins), np.empty(Nbins), np.empty(Nbins), np.empty(Nbins)
    for i in range(Nbins):
        x[i], y[i], x_err[i], y_err[i]  = hist.GetBinCenter(i+1), hist.GetBinContent(i+1), hist.GetBinWidth(i+1), hist.GetBinError(i+1)
    return [x, y, x_err, y_err]

def fitfunc(x, c1, sigma1, mu, a, b):
    return c1/sigma1/np.sqrt(2*np.pi)*np.exp(-(x-mu)**2/(2*sigma1**2)) + a*x + b


edges_energy = np.arange(6.0, 11.0, 0.5)        
edges_theta = np.arange(20, 165, 5)
edges_theta

file_data   = root.TFile("../data.root")
tree_data   = file_data.Get("flattree_piminus_p_4He_recon")
yield_hist_data = create_yield_hist()

for i in range(tree_data.GetEntries()):
    tree_data.GetEntry(i)
    MissingPMinus = tree_data.MissingPMinus
    BeamEnergy    = tree_data.BeamEnergy
    thetaCM       = tree_data.thetaCM
    WeightFactor  = tree_data.WeightFactor
    
    for j in range(len(edges_energy)-1):
        if BeamEnergy < edges_energy[j+1]:
            for k in range(len(edges_theta)-1):
                if thetaCM < edges_theta[k+1]:
                    yield_hist_data[j][k].Fill(MissingPMinus, WeightFactor)
                    break
            break

yield_data = np.zeros((len(edges_energy)-1, len(edges_theta)-1))
error_data = np.zeros((len(edges_energy)-1, len(edges_theta)-1))

with PdfPages('output/fit_piminus_p_4He_data.pdf') as pdf:

    for i in range(len(edges_energy)-1):        
        for j in range(len(edges_theta)-1):
            
#             points = hist_points_1D(yield_hist_data[i][j])

#             popt, pcov = curve_fit(gaussian, points[0][points[1].nonzero()], points[1][points[1].nonzero()], sigma=points[3][points[1].nonzero()], absolute_sigma=True, p0=[points[1].sum()/100, 0.4, 0.95], bounds=[[0, 0.1, 0.5], [points[1].sum(), 1, 1.3]])
#             perr = np.sqrt(np.diag(pcov))
            
#             yield_data_energy.append(popt[0]/hist_norm)
#             yield_data_err_energy.append(perr[0]/hist_norm)

#             fig, ax = plt.subplots(figsize=(10, 7))
#             plt.errorbar(points[0], points[1], xerr=points[2], yerr=points[3], fmt='o')
#             plt.plot(np.linspace(-2, 4, 1000), gaussian(np.linspace(-2, 4, 1000), popt[0], popt[1], popt[2]))
#             plt.xlabel('$M^2_{miss} [GeV^2]$')
#             plt.ylabel('Counts/0.02 $GeV^2$')
#             plt.text(1, 1, \
#                      '${0:.3f} \, GeV \, < \, E_\\gamma \, < \, {1:.3f} \, GeV$ \n \
#                       ${2:.1f} \, GeV^2 \, < \, -t \, < \, {3:.1f} \, GeV^2$ \n \
#                       $N_{{\\pi^-}} \, = \, {4:.3f} \, \pm \, {5:.3f}$ \n \
#                       $\sigma \, = \, {6:.3f} \, \pm \, {7:.3f} \, GeV^2$ \n \
#                       $M^2_p \, = \, {8:.3f} \, \pm \, {9:.3f} \, GeV^2$'\
#                      .format(edges_energy[i], edges_energy[i+1], edges_minust[i][j], edges_minust[i][j+1], popt[0]/hist_norm, perr[0]/hist_norm, popt[1], perr[1], popt[2], perr[2]), \
#                      fontsize=15, ha='right', va='top', transform=ax.transAxes\
#                     )

#             pdf.savefig()
#             plt.close()
        
#         yield_data.append(yield_data_energy)
#         yield_data_err.append(yield_data_err_energy)

points = hist_points_1D(yield_hist_data[3][3])
plt.errorbar(points[0], points[1], xerr=points[2], yerr=points[3], fmt='o')

popt, pcov = curve_fit(fitfunc, points[0][points[1].nonzero()], points[1][points[1].nonzero()], sigma=points[3][points[1].nonzero()], absolute_sigma=True, p0=[points[1].sum(), 0.4, 0.95, -200, 200] , bounds=[[100, 0.1, 0.5, -1000, 0], [points[1].sum(), 1, 1.3, 0, 1000]])
# err = np.sqrt(np.diag(pcov))     
plt.plot(np.linspace(0.5, 1.3, 1000), fitfunc(np.linspace(0.5, 1.3, 1000), popt[0], popt[1], popt[2], popt[3], popt[4]))

file_sim    = root.TFile("/work/halld2/home/boyu/src_analysis/selection/output/flattree_piminus_p_2H_sim.root")
tree_sim    = file_sim.Get("piminus_p_2H_recon")
yield_hist_sim = create_yield_hist(edges_minust)

for i in range(tree_sim.GetEntries()):
    tree_sim.GetEntry(i)
    BeamP4_Measured    = tree_sim.BeamP4_Measured
    PiMinusP4_Measured = tree_sim.PiMinusP4_Measured
    ProtonP4_Measured  = tree_sim.ProtonP4_Measured
    MissingP4_Measured = tree_sim.MissingP4_Measured
    AccidWeight      = tree_sim.HistAccidWeightFactor
    L1TriggerBits    = tree_sim.L1TriggerBits
    
    MinusT_Measured = -(BeamP4_Measured - PiMinusP4_Measured).Mag2()
    
    if L1TriggerBits == 0:
        continue
    
    for j in range(len(edges_energy)-1):
        if BeamP4_Measured.E() < edges_energy[j+1]:
            for k in range(len(edges_minust[j])-1):
                if MinusT_Measured < edges_minust[j][k+1]:
                    yield_hist_sim[j][k].Fill(MissingP4_Measured.M2(), AccidWeight)
                    break
            break

yield_sim = []
yield_sim_err = []

from matplotlib.backends.backend_pdf import PdfPages
with PdfPages('output/fit_piminus_p_2H_sim.pdf') as pdf:

    for i in range(len(edges_energy)-1):
        
        yield_sim_energy = []
        yield_sim_err_energy = []
        
        for j in range(len(edges_minust[i])-1):
            
            points = hist_points_1D(yield_hist_sim[i][j])

            popt, pcov = curve_fit(gaussian, points[0][points[1].nonzero()], points[1][points[1].nonzero()], sigma=points[3][points[1].nonzero()], absolute_sigma=True, p0=[points[1].sum()/100, 0.4, 0.95], bounds=[[0, 0.1, 0.5], [points[1].sum(), 1, 1.3]])
            perr = np.sqrt(np.diag(pcov))
            
            yield_sim_energy.append(popt[0]/hist_norm)
            yield_sim_err_energy.append(perr[0]/hist_norm)

            fig, ax = plt.subplots(figsize=(10, 7))
            plt.errorbar(points[0], points[1], xerr=points[2], yerr=points[3], fmt='o')
            plt.plot(np.linspace(-2, 4, 1000), gaussian(np.linspace(-2, 4, 1000), popt[0], popt[1], popt[2]))
            plt.xlabel('$M^2_{miss} [GeV^2]$')
            plt.ylabel('Counts/0.02 $GeV^2$')
            plt.text(1, 1, \
                     '${0:.3f} \, GeV \, < \, E_\\gamma \, < \, {1:.3f} \, GeV$ \n \
                      ${2:.1f} \, GeV^2 \, < \, -t \, < \, {3:.1f} \, GeV^2$ \n \
                      $N_{{\\pi^-}} \, = \, {4:.3f} \, \pm \, {5:.3f}$ \n \
                      $\sigma \, = \, {6:.3f} \, \pm \, {7:.3f} \, GeV^2$ \n \
                      $M^2_p \, = \, {8:.3f} \, \pm \, {9:.3f} \, GeV^2$'\
                     .format(edges_energy[i], edges_energy[i+1], edges_minust[i][j], edges_minust[i][j+1], popt[0]/hist_norm, perr[0]/hist_norm, popt[1], perr[1], popt[2], perr[2]), \
                     fontsize=15, ha='right', va='top', transform=ax.transAxes\
                    )

            pdf.savefig()
            plt.close()
            
        yield_sim.append(yield_sim_energy)
        yield_sim_err.append(yield_sim_err_energy)

file_thrown = root.TFile("/work/halld2/home/boyu/src_analysis/selection/output/flattree_piminus_p_2H_thrown.root")
tree_thrown = file_thrown.Get("piminus_p_2H_thrown")
yield_thrown, yield_thrown_err = create_zero_list(edges_minust), create_zero_list(edges_minust)

for i in range(tree_thrown.GetEntries()):
    tree_thrown.GetEntry(i)
    BeamP4_Thrown    = tree_thrown.BeamP4_Thrown
    PiMinusP4_Thrown = tree_thrown.PiMinusP4_Thrown
    ProtonP4_Thrown  = tree_thrown.ProtonP4_Thrown
    MissingP4_Thrown = tree_thrown.MissingP4_Thrown
    
    MinusT_Thrown = -(BeamP4_Thrown - PiMinusP4_Thrown).Mag2()
    
    for j in range(len(edges_energy)-1):
        if BeamP4_Thrown.E() < edges_energy[j+1]:
            for k in range(len(edges_minust[j])-1):
                if MinusT_Thrown < edges_minust[j][k+1]:
                    yield_thrown[j][k] += 1
                    yield_thrown_err[j][k] = np.sqrt(yield_thrown[j][k])
                    break
            break

acceptance, acceptance_error = create_zero_list(edges_minust), create_zero_list(edges_minust)
for i in range(len(edges_energy)-1):
    for j in range(len(edges_minust[i])-1):
        acceptance[i][j] = yield_sim[i][j]/yield_thrown[i][j]
        acceptance_error[i][j] = acceptance[i][j]*(yield_sim_err[i][j]/yield_sim[i][j] + yield_thrown_err[i][j]/yield_thrown[i][j])

lumi = [18.5985329-16.4927209, 16.4927209-15.2246128, 15.2246128-13.5462486, 13.5462486-10.2872631, 10.2872631-3.9112386, 3.9112386]

cs, cs_error = create_zero_list(edges_minust), create_zero_list(edges_minust)
for i in range(len(edges_energy)-1):
    for j in range(len(edges_minust[i])-1):
        cs[i][j] = yield_data[i][j]/acceptance[i][j]/(edges_minust[i][j+1] - edges_minust[i][j])/lumi[i]/1000
        cs_error[i][j] = cs[i][j]*(yield_data_err[i][j]/yield_data[i][j] + acceptance_error[i][j]/acceptance[i][j])


#Constants
mass_neutron = 0.93956542052
mass_proton  = 0.93827208816
mass_piminus = 0.13957039
mass_photon  = 0
rad_to_deg = 180/3.1415926535897932

#Function to convert the square root of s to the photon energy, assuming neutron at rest
def SqrtS_to_Eg(SqrtS):
    return (SqrtS**2-mass_neutron**2)/(2*mass_neutron)

#Function to convert the photon energy to the square root of s, assuming neutron at rest
def Eg_to_SqrtS(Eg):
    return np.sqrt(mass_neutron**2+2*mass_neutron*Eg)

#Function to calculate the upper limit of minus t, given a specific center of mass energy squared s
def abst_max(s):
    return abs(mass_photon**2+mass_neutron**2+mass_piminus**2+mass_proton**2-s)

#Function to calculate the charged pion ratio, given s and t
def piminus_to_piplus_ratio(s,t):
    u = mass_photon**2+mass_neutron**2+mass_piminus**2+mass_proton**2 - s - t
    e_d = -1/3
    e_u = 2/3
    return ((e_d*u+e_u*s)/(e_u*u+e_d*s))**2

plt.rcParams['figure.dpi'] = 200
plt.xlabel("$-t/GeV^2$")
plt.ylabel("$d\sigma/dt (nb/GeV^2)$")
plt.yscale('log')
plt.title("$E_\gamma$ = 7.5 GeV")
plt.ylim(1e-2, 1e3)
plt.xlim(0,14)

points1 = []
points2 = []
for i in range(len(edges_minust[2])-1):
    points1.append((edges_minust[2][i]+edges_minust[2][i+1])/2)
    points2.append((edges_minust[2][i]-edges_minust[2][i+1])/2)

plt.errorbar(points1[:-1], cs[2][:-1], xerr = points2[:-1], yerr = cs_error[2][:-1], marker='.', fmt='o', label="$\pi^- p $ (This work)")
plt.yscale('log')

points_anderson = np.array([1.95, 2.24, 2.54, 3.02, 4.12, 5.04, 6.17, 6.61, 7.17, 8.23, 9.40, 10.43, 11.26, 11.63])
cs_anderson = np.array([3.45, 1.15, 0.481, 0.212, 0.121, 0.065, 0.050, 0.045, 0.03, 0.059, 0.03, 0.169, 0.374, 0.83])
cs_error_anderson = np.array([0.16, 0.12, 0.079, 0.051, 0.019, 0.023, 0.013, 0.019, 0.018, 0.033, 0.026, 0.038, 0.097, 0.14])
pion_ratio = piminus_to_piplus_ratio(np.ones(14)*Eg_to_SqrtS(7.5)**2, -points_anderson)
plt.errorbar(points_anderson, cs_anderson*pion_ratio, yerr=cs_error_anderson*pion_ratio, marker='.', fmt='o', label="scaled $\pi^+ n $ (Anderson)")

plt.legend()